# -*- coding: utf-8 -*-
"""
engines/toolchain.py -- resolucao de EXECUTAVEIS externos (Vina, OpenBabel,
Boltz, HADDOCK3, nvidia-smi).

Distinto de _deps.py (que cuida de pacotes Python via pip). Aqui o problema e o
acoplamento historico a caminhos absolutos espalhados pelo codigo: vina.exe em
Projetos\\docking\\, obabel no env 'datascience', boltz no env 'foldgpu'.
`find_tool()` centraliza isso em UMA cascata previsivel:

    1. PATH               (shutil.which)
    2. registry/toolchain.json     (override declarativo, versionavel)
    3. heuristica conda   (<conda_root>/envs/<env>/{,'Scripts','Library/bin','bin'})

Se nada resolver, devolve None -- e o chamador (engine.check_toolchain / run)
falha com mensagem clara. NUNCA adivinha um binario plausivel nem instala nada.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _toolchain_overrides():
    """Le registry/toolchain.json, se existir: {tool: caminho_ou_{path:...}}."""
    try:
        import registry
        path = registry.REGISTRY_DIR / "toolchain.json"
    except Exception:
        path = ROOT / "registry" / "toolchain.json"
    if not Path(path).exists():
        return {}
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return {}


def _conda_root():
    """Raiz da instalacao conda, inferida do interpretador atual."""
    exe = Path(sys.executable).resolve()
    parts = exe.parts
    if "envs" in parts:
        return Path(*parts[: parts.index("envs")])
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        p = Path(prefix)
        return p.parent.parent if p.parent.name == "envs" else p
    # fallback historico do projeto
    guess = Path.home() / "miniforge3"
    return guess if guess.exists() else exe.parent


def _exe_candidates(name):
    """Nomes de arquivo possiveis para um executavel (com/sem .exe no Windows)."""
    names = [name]
    if os.name == "nt" and not name.lower().endswith(".exe"):
        names.append(name + ".exe")
    return names


def _search_conda(exe_names, conda_envs):
    root = _conda_root()
    subdirs = ("", "Scripts", os.path.join("Library", "bin"), "bin")
    for env in conda_envs:
        env_dir = root / "envs" / env
        for sub in subdirs:
            for exe in exe_names:
                cand = env_dir / sub / exe if sub else env_dir / exe
                if cand.exists():
                    return str(cand)
    return None


def find_tool(name, *, exe_names=None, conda_envs=(), registry_key=None,
              extra_paths=()):
    """Resolve o caminho de um executavel pela cascata PATH -> registry -> conda.

    name         : nome logico (ex.: 'vina', 'obabel', 'boltz', 'haddock3').
    exe_names    : nomes de arquivo alternativos (default: [name] + name.exe).
    conda_envs   : envs conda onde procurar (ex.: ['datascience'], ['foldgpu']).
    registry_key : chave em registry/toolchain.json (default: name).
    extra_paths  : caminhos absolutos extras a testar por ultimo.
    Devolve o caminho (str) ou None. NUNCA levanta.
    """
    exe_names = exe_names or _exe_candidates(name)

    # 1. PATH
    for exe in exe_names:
        hit = shutil.which(exe)
        if hit:
            return hit

    # 2. override declarativo em registry/toolchain.json
    overrides = _toolchain_overrides()
    entry = overrides.get(registry_key or name)
    if isinstance(entry, dict):
        entry = entry.get("path")
    if entry and Path(entry).exists():
        return str(entry)

    # 3. heuristica conda
    hit = _search_conda(exe_names, conda_envs)
    if hit:
        return hit

    # 4. caminhos extras explicitos
    for p in extra_paths:
        if p and Path(p).exists():
            return str(p)

    return None


def sha256_file(path):
    """sha256 hex de um arquivo (blocos). None se ilegivel. Nunca levanta."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def tool_version(path, args=("--version",), timeout=20):
    """Primeira linha da saida de `path --version` (ou stderr). None se falhar."""
    try:
        r = subprocess.run([str(path), *args], capture_output=True, text=True,
                           timeout=timeout)
        blob = (r.stdout or r.stderr or "").strip()
        return blob.splitlines()[0][:200] if blob else None
    except Exception:
        return None


def tool_entry(path, *, version_args=("--version",), hash_binary=True):
    """{available, path, sha256?, size_bytes?, version} de um executavel."""
    if not path:
        return {"available": False}
    p = Path(path)
    entry = {"available": True, "path": str(p)}
    if p.exists():
        if hash_binary:
            entry["sha256"] = sha256_file(p)
        try:
            entry["size_bytes"] = p.stat().st_size
        except OSError:
            pass
    entry["version"] = tool_version(p, version_args)
    return entry


def gpu_entry():
    """Estado da GPU via nvidia-smi. {available, name, driver_version, memory_total}."""
    try:
        r = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,driver_version,memory.total",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=20)
        rows = (r.stdout or "").strip().splitlines()
        if rows:
            name, driver, mem = (x.strip() for x in rows[0].split(",")[:3])
            return {"available": True, "name": name,
                    "driver_version": driver, "memory_total": mem}
    except Exception as exc:
        return {"available": False, "error": repr(exc)[:200]}
    return {"available": False}


# Especificacao dos executaveis conhecidos da plataforma: onde procurar cada um.
# Alimenta find_tool() e o snapshot de toolchain do runner.
KNOWN_TOOLS = {
    "vina":     {"conda_envs": (), "exe_names": None},
    "obabel":   {"conda_envs": ("datascience",), "exe_names": None},
    "boltz":    {"conda_envs": ("foldgpu",), "exe_names": None},
    "haddock3": {"conda_envs": ("haddock3", "datascience"), "exe_names": None},
}


def resolve_known(name):
    """find_tool para um tool conhecido, usando KNOWN_TOOLS. None se desconhecido."""
    spec = KNOWN_TOOLS.get(name)
    if spec is None:
        return find_tool(name)
    return find_tool(name, exe_names=spec.get("exe_names"),
                     conda_envs=spec.get("conda_envs", ()))


_SNAPSHOT_CACHE = {}


def snapshot(tools=("vina", "obabel", "boltz", "haddock3"), force=False):
    """Snapshot da toolchain para proveniencia: versao+hash de cada tool + GPU.

    Cacheado em memoria por processo: sondar os binarios (subprocess `--version`,
    hash) custa ~10 s -- `boltz --version` carrega o pacote inteiro -- e nao muda
    durante a vida do servidor. O runner embute este snapshot em CADA run.json;
    sem cache, todo run (inclusive dry_run) pagaria os 10 s de novo, e o run.json
    so nasceria depois disso. Passe force=True para re-sondar (ex.: apos instalar
    um binario que faltava)."""
    key = tuple(tools)
    if not force and key in _SNAPSHOT_CACHE:
        return _SNAPSHOT_CACHE[key]
    out = {"python": {"available": True, "version": sys.version.split()[0],
                      "executable": sys.executable}}
    for name in tools:
        out[name] = tool_entry(resolve_known(name))
    out["gpu"] = gpu_entry()
    _SNAPSHOT_CACHE[key] = out
    return out
