# -*- coding: utf-8 -*-
"""tools/deploy_hostinger.py -- deploy do viewer estatico para o Hostinger por
FTP(S)/SFTP, sem clicar no hPanel.

SEGURANCA: as credenciais NUNCA vivem no codigo nem no git. Ficam num arquivo
LOCAL gitignored `.hostinger.json` (ou em variaveis de ambiente HOSTINGER_*).
Quem preenche a senha e o AUTOR, no arquivo dele; este script so a LE.

.hostinger.json (na raiz do repo, gitignored) -- modelo em .hostinger.example.json:
{
  "protocol": "ftps",              // "ftps" (padrao, Hostinger shared) ou "sftp"
  "host": "ftp.SEUDOMINIO",        // hPanel -> Files -> FTP Accounts (host/usuario)
  "port": 21,                       // 21 p/ ftps, 22 p/ sftp
  "user": "u123456789",
  "password": "SUA_SENHA_FTP",
  "remote_dir": "public_html/docking"
}

Uso (DRY-RUN por padrao -- nada e enviado/apagado sem --apply):
  python tools/deploy_hostinger.py --list
  python tools/deploy_hostinger.py --push --local <destino/docking_suite_lgm> --apply
  python tools/deploy_hostinger.py --clean-orphans --apply
  python tools/deploy_hostinger.py --push --local <pkg> --clean-orphans --apply   # tudo
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Orfaos do deploy antigo (single-campanha) a remover em public_html/docking/.
# A landing nova (index.html) + SUB01/ + SUB04/ nao referenciam nada disto.
ORPHAN_FILES = ["config.json", "CHANGELOG.md", "README_DEPLOY.md", "README_HOSTINGER.md"]
ORPHAN_DIRS = ["data", "figures", "assets", "docs", "lib"]
# NUNCA remover: index.html, SUB01/, SUB04/, _auth_path.php, .htaccess


def load_cfg():
    """Le .hostinger.json (gitignored) e/ou env HOSTINGER_*. Env sobrepoe arquivo."""
    cfg = {}
    p = ROOT / ".hostinger.json"
    if p.exists():
        cfg.update(json.loads(p.read_text(encoding="utf-8")))
    for k in ("protocol", "host", "port", "user", "password", "remote_dir"):
        v = os.environ.get("HOSTINGER_" + k.upper())
        if v:
            cfg[k] = v
    cfg.setdefault("protocol", "ftps")
    cfg.setdefault("port", 21 if cfg["protocol"] == "ftps" else 22)
    cfg.setdefault("remote_dir", "public_html/docking")
    # Sanitiza o host: aceita colado com esquema/barra (ftp://host/, sftp://host).
    h = str(cfg.get("host", "")).strip()
    for pre in ("ftps://", "ftp://", "sftp://", "ssh://", "https://", "http://"):
        if h.lower().startswith(pre):
            h = h[len(pre):]
    cfg["host"] = h.strip("/").split("/")[0]
    missing = [k for k in ("host", "user", "password") if not cfg.get(k)]
    if missing:
        sys.exit(f"[erro] credenciais faltando: {missing}. Crie .hostinger.json "
                 f"(veja .hostinger.example.json) ou defina HOSTINGER_HOST/USER/PASSWORD.")
    return cfg


# ---------------------------------------------------------------- backend FTPS
class FTPSBackend:
    def __init__(self, cfg):
        from ftplib import FTP_TLS
        self.ftp = FTP_TLS()
        self.ftp.connect(cfg["host"], int(cfg["port"]), timeout=30)
        self.ftp.login(cfg["user"], cfg["password"])
        self.ftp.prot_p()               # protege o canal de dados

    def listdir(self, path):
        try:
            return self.ftp.nlst(path)
        except Exception:
            return []

    def is_dir(self, path):
        cur = self.ftp.pwd()
        try:
            self.ftp.cwd(path); self.ftp.cwd(cur); return True
        except Exception:
            return False

    def mkdirs(self, path):
        parts = path.strip("/").split("/")
        cur = ""
        for p in parts:
            cur += "/" + p
            try:
                self.ftp.mkd(cur)
            except Exception:
                pass

    def upload(self, local: Path, remote: str):
        with open(local, "rb") as fh:
            self.ftp.storbinary("STOR " + remote, fh)

    def delete_file(self, remote: str):
        self.ftp.delete(remote)

    def rmdir_recursive(self, remote: str):
        for name in list(self.listdir(remote)):
            base = name.split("/")[-1]
            if base in (".", ".."):
                continue
            full = remote + "/" + base
            if self.is_dir(full):
                self.rmdir_recursive(full)
            else:
                try: self.ftp.delete(full)
                except Exception: pass
        try: self.ftp.rmd(remote)
        except Exception: pass

    def close(self):
        try: self.ftp.quit()
        except Exception: pass


# ---------------------------------------------------------------- backend SFTP
class SFTPBackend:
    def __init__(self, cfg):
        import paramiko
        self.t = paramiko.Transport((cfg["host"], int(cfg["port"])))
        self.t.connect(username=cfg["user"], password=cfg["password"])
        self.sftp = paramiko.SFTPClient.from_transport(self.t)

    def listdir(self, path):
        try: return self.sftp.listdir(path)
        except Exception: return []

    def is_dir(self, path):
        import stat
        try: return stat.S_ISDIR(self.sftp.stat(path).st_mode)
        except Exception: return False

    def mkdirs(self, path):
        cur = ""
        for p in path.strip("/").split("/"):
            cur += "/" + p
            try: self.sftp.mkdir(cur)
            except Exception: pass

    def upload(self, local: Path, remote: str):
        self.sftp.put(str(local), remote)

    def delete_file(self, remote: str):
        self.sftp.remove(remote)

    def rmdir_recursive(self, remote: str):
        for base in list(self.listdir(remote)):
            full = remote + "/" + base
            if self.is_dir(full): self.rmdir_recursive(full)
            else:
                try: self.sftp.remove(full)
                except Exception: pass
        try: self.sftp.rmdir(remote)
        except Exception: pass

    def close(self):
        try: self.sftp.close(); self.t.close()
        except Exception: pass


def backend(cfg):
    return SFTPBackend(cfg) if cfg["protocol"] == "sftp" else FTPSBackend(cfg)


def push_tree(be, local_dir: Path, remote_dir: str, apply: bool):
    n = 0
    for f in sorted(local_dir.rglob("*")):
        if f.is_dir():
            continue
        rel = f.relative_to(local_dir).as_posix()
        remote = remote_dir + "/" + rel
        print(f"  UPLOAD {rel}")
        if apply:
            be.mkdirs(remote.rsplit("/", 1)[0])
            be.upload(f, remote)
        n += 1
    print(f"[{'ENVIADO' if apply else 'dry-run'}] {n} arquivos -> {remote_dir}")


def clean_orphans(be, remote_dir: str, apply: bool):
    for d in ORPHAN_DIRS:
        full = remote_dir + "/" + d
        if be.is_dir(full):
            print(f"  RM DIR  {d}/")
            if apply: be.rmdir_recursive(full)
    for f in ORPHAN_FILES:
        full = remote_dir + "/" + f
        if be.listdir(full) or not be.is_dir(full):
            # tenta apagar; se nao existir, ignora
            print(f"  RM FILE {f}")
            if apply:
                try: be.delete_file(full)
                except Exception: pass
    print(f"[{'LIMPO' if apply else 'dry-run'}] orfaos em {remote_dir}")


def main():
    ap = argparse.ArgumentParser(description="Deploy do viewer estatico p/ Hostinger (FTP/SFTP).")
    ap.add_argument("--local", default=str(Path(os.environ.get(
        "DOCKING_LGM_DIR",
        r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_suite_lgm"))),
        help="pasta local do pacote (docking_suite_lgm)")
    ap.add_argument("--push", action="store_true", help="envia o pacote local")
    ap.add_argument("--clean-orphans", action="store_true", help="remove os orfaos do deploy antigo")
    ap.add_argument("--list", action="store_true", help="lista o diretorio remoto")
    ap.add_argument("--apply", action="store_true", help="EXECUTA (sem isso e dry-run)")
    args = ap.parse_args()

    cfg = load_cfg()
    rd = cfg["remote_dir"]
    print(f"[cfg] {cfg['protocol']}://{cfg['user']}@{cfg['host']}:{cfg['port']}/{rd}  "
          f"({'APPLY' if args.apply else 'DRY-RUN'})")
    be = backend(cfg)
    try:
        if args.list or not (args.push or args.clean_orphans):
            print("[remoto]", sorted(x.split('/')[-1] for x in be.listdir(rd)))
        if args.clean_orphans:
            clean_orphans(be, rd, args.apply)
        if args.push:
            local = Path(args.local)
            if not local.is_dir():
                sys.exit(f"[erro] pasta local nao existe: {local}")
            push_tree(be, local, rd, args.apply)
    finally:
        be.close()
    print("OK." + ("" if args.apply else "  (dry-run -- rode com --apply p/ efetivar)"))


if __name__ == "__main__":
    main()
