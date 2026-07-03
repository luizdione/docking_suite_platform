#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""enable_basic_auth.py - ativa protecao por senha (Basic Auth) num bundle dist_web/.

Gera `<out>/.htpasswd` (apenas o HASH bcrypt da senha; a senha em texto nunca e
gravada) e reescreve `<out>/.htaccess` ativando o bloco Basic Auth. Cria tambem um
helper `_auth_path.php` que, acessado uma vez no servidor, imprime o caminho absoluto
correto para o AuthUserFile.

Esta ferramenta NAO contem credenciais: usuario e senha vem por argumento e o
`.htpasswd`/senha NUNCA sao versionados (o dist_web/ e gitignored).

Uso:
    python enable_basic_auth.py --out dist_web --user darwin --password "SENHA" \
        --htpasswd-path "/home/uXXXX/domains/<dominio>/public_html/docking/.htpasswd"

Se --htpasswd-path for omitido, usa um caminho provavel e deixa instrucoes no .htaccess.
Alternativa sem arquivos: hPanel -> Avancado -> "Protecao por senha de diretorios".
"""
from __future__ import annotations

import argparse
from pathlib import Path

import bcrypt

DEFAULT_PATH = ("/home/u112348440/domains/ifrj-crj-geneticamolecular.online/"
                "public_html/docking/.htpasswd")

HTACCESS_BASE = """# =============================================================================
# NQ-Docking viewer - Apache/LiteSpeed (Hostinger). Basic Auth ATIVADO por
# enable_basic_auth.py. O .htpasswd guarda apenas o hash bcrypt da senha.
# =============================================================================
AddDefaultCharset UTF-8
AddType application/json .json
AddType image/svg+xml .svg
AddType text/markdown .md
AddType application/wasm .wasm
AddType application/pdf .pdf
Options -Indexes

<IfModule mod_expires.c>
  ExpiresActive On
  ExpiresByType application/javascript "access plus 7 days"
  ExpiresByType image/png "access plus 7 days"
</IfModule>

# --- Basic Auth -------------------------------------------------------------
AuthType Basic
AuthName "NQ-Docking - area restrita (IFRJ)"
AuthUserFile "{HTPASSWD_PATH}"
Require valid-user

# Nao permitir baixar os proprios arquivos de credencial/segredo.
<FilesMatch "^\\.(htpasswd|htaccess)$">
  Require all denied
</FilesMatch>
"""

AUTH_PATH_PHP = """<?php
// Helper temporario: acesse esta pagina UMA vez para descobrir o caminho absoluto
// desta pasta no servidor, cole-o no AuthUserFile do .htaccess e APAGUE este arquivo.
header('Content-Type: text/plain; charset=utf-8');
echo "Diretorio absoluto desta pasta:\\n";
echo __DIR__ . "\\n\\n";
echo "AuthUserFile para o .htaccess:\\n";
echo '  AuthUserFile "' . __DIR__ . '/.htpasswd"' . "\\n";
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="dist_web")
    ap.add_argument("--user", required=True)
    ap.add_argument("--password", required=True)
    ap.add_argument("--htpasswd-path", default=DEFAULT_PATH,
                    help="caminho ABSOLUTO do .htpasswd no servidor (para o AuthUserFile)")
    args = ap.parse_args()

    out = Path(args.out)
    if not out.exists():
        raise SystemExit(f"bundle nao encontrado: {out} (rode export_web.py antes)")

    # 1) .htpasswd (apenas hash bcrypt; formato apache: user:$2y$...)
    salted = bcrypt.hashpw(args.password.encode("utf-8"), bcrypt.gensalt(rounds=10))
    # Apache/LiteSpeed reconhecem o prefixo $2y$; o bcrypt do Python emite $2b$ (equivalente).
    hash_str = salted.decode("ascii").replace("$2b$", "$2y$", 1)
    (out / ".htpasswd").write_text(f"{args.user}:{hash_str}\n", encoding="ascii")

    # 2) .htaccess ativado
    (out / ".htaccess").write_text(
        HTACCESS_BASE.replace("{HTPASSWD_PATH}", args.htpasswd_path), encoding="utf-8")

    # 3) helper para descobrir o caminho absoluto
    (out / "_auth_path.php").write_text(AUTH_PATH_PHP, encoding="utf-8")

    print(f"[OK] Basic Auth ativado em {out}/")
    print(f"     usuario: {args.user}  (hash bcrypt gravado; senha NAO gravada)")
    print(f"     AuthUserFile: {args.htpasswd_path}")
    print("     Se der erro 401/500: acesse _auth_path.php no site para pegar o caminho")
    print("     exato, ajuste o AuthUserFile no .htaccess e apague _auth_path.php.")


if __name__ == "__main__":
    main()
