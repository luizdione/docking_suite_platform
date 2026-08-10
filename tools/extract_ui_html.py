# -*- coding: utf-8 -*-
"""
Extrai a UI de dashboard_v2.py para static/app.html.

A interface inteira vivia como uma string Python de ~275 mil caracteres dentro
do servidor. Isso tornava qualquer mudanca visual um diff num arquivo de 374 KB,
e impedia editar HTML/CSS/JS com as ferramentas de HTML/CSS/JS.

A extracao e MECANICA: o conteudo entre os delimitadores da string vai para o
arquivo byte a byte, e o servidor passa a le-lo. Nenhum caractere e alterado.
O gate e o SHA-256 do HTML servido, comparado antes e depois.

Jinja continua fora de questao: `{{ }}` colide com os template literals do JS
embutido. O arquivo e servido cru, com a mesma substituicao de placeholders.

Uso:
    python tools/extract_ui_html.py --check     # so verifica, nao escreve
    python tools/extract_ui_html.py
"""
import argparse
import hashlib
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DASHBOARD = ROOT / "dashboard_v2.py"
STATIC_DIR = ROOT / "static"
APP_HTML = STATIC_DIR / "app.html"

REPLACEMENT = '''# ─── UI ───────────────────────────────────────────────────────────
# A interface vive em static/app.html (extraida de uma string de ~275 mil
# caracteres que ficava aqui). E servida CRUA, sem motor de template: Jinja
# colidiria com os template literals `${...}` do JS embutido. A injecao de
# dados continua sendo substituicao de placeholder, como antes.
_STATIC_DIR = Path(__file__).resolve().parent / "static"
_APP_HTML = _STATIC_DIR / "app.html"
_HTML_CACHE = {"mtime": None, "text": None}


def _load_html():
    """Le a UI do disco, com cache invalidado por mtime.

    Recarregar em desenvolvimento sem reiniciar o servidor e o ganho pratico
    imediato da extracao.
    """
    try:
        mtime = _APP_HTML.stat().st_mtime_ns
    except OSError as exc:
        raise RuntimeError(
            f"UI nao encontrada em {_APP_HTML}. O arquivo static/app.html faz parte "
            "do repositorio; se sumiu, restaure-o do git."
        ) from exc
    if _HTML_CACHE["mtime"] != mtime:
        _HTML_CACHE["text"] = _APP_HTML.read_text(encoding="utf-8")
        _HTML_CACHE["mtime"] = mtime
    return _HTML_CACHE["text"]


class _HTMLProxy(str):
    """Compatibilidade: o codigo existente usa `HTML.replace(...)`.

    Mantem `HTML` como se fosse a string de antes, mas relendo do disco a cada
    acesso. Assim nenhum chamador precisa mudar.
    """

    def __new__(cls):
        return super().__new__(cls, "")

    def replace(self, *args, **kwargs):
        return _load_html().replace(*args, **kwargs)

    def __str__(self):
        return _load_html()

    def __len__(self):
        return len(_load_html())

    def __contains__(self, item):
        return item in _load_html()


HTML = _HTMLProxy()
'''


def sha256(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def find_html_literal(src):
    """Localiza a string HTML. Devolve (inicio_stmt, inicio_conteudo, fim_conteudo, fim_stmt)."""
    match = re.search(r"^HTML\s*=\s*(r?)(\"\"\"|''')", src, re.M)
    if not match:
        return None
    delim = match.group(2)
    content_start = match.end()
    content_end = src.index(delim, content_start)
    return match.start(), content_start, content_end, content_end + len(delim)


def current_served_html():
    """SHA-256 do HTML que o servidor entrega hoje em '/'."""
    sys.path.insert(0, str(ROOT))
    import dashboard_v2

    client = dashboard_v2.app.test_client()
    return client.get("/").data.decode("utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="nao escreve, so relata")
    args = ap.parse_args()

    src = DASHBOARD.read_text(encoding="utf-8")
    found = find_html_literal(src)
    if not found:
        if APP_HTML.exists():
            print("[ok] extracao ja aplicada: dashboard_v2.py nao tem mais a string HTML")
            return 0
        print("[erro] nao encontrei a string HTML em dashboard_v2.py")
        return 1

    stmt_start, content_start, content_end, stmt_end = found
    html_text = src[content_start:content_end]
    print(f"string HTML: {len(html_text)} chars "
          f"(linhas {src[:stmt_start].count(chr(10))+1}-{src[:content_end].count(chr(10))+1})")

    print("capturando o HTML servido ANTES da mudanca...")
    before = current_served_html()
    print(f"  sha256 antes: {sha256(before)}  ({len(before)} chars)")

    if args.check:
        print("(--check: nada escrito)")
        return 0

    backup = DASHBOARD.with_suffix(".py.bak_ui_extract")
    shutil.copyfile(DASHBOARD, backup)
    print(f"backup: {backup.name}")

    STATIC_DIR.mkdir(exist_ok=True)
    # newline='' para nao traduzir as quebras de linha: qualquer normalizacao
    # aqui muda o hash e, com ele, a prova de que a extracao foi neutra.
    with open(APP_HTML, "w", encoding="utf-8", newline="") as fh:
        fh.write(html_text)
    print(f"escrito: static/app.html ({len(html_text)} chars)")

    new_src = src[:stmt_start] + REPLACEMENT + src[stmt_end:]
    with open(DASHBOARD, "w", encoding="utf-8", newline="") as fh:
        fh.write(new_src)
    print(f"dashboard_v2.py: {len(src)} -> {len(new_src)} chars "
          f"({(len(src)-len(new_src))/1024:.0f} KB a menos)")

    print("\nrecarregando e comparando...")
    for mod in ("dashboard_v2",):
        sys.modules.pop(mod, None)
    after = current_served_html()
    h_before, h_after = sha256(before), sha256(after)
    print(f"  sha256 antes : {h_before}")
    print(f"  sha256 depois: {h_after}")

    if h_before == h_after:
        print("\n[OK] HTML servido identico byte a byte. Extracao neutra.")
        return 0

    print("\n[FALHA] o HTML servido MUDOU. Revertendo.")
    shutil.copyfile(backup, DASHBOARD)
    print(f"  dashboard_v2.py restaurado do backup")
    print(f"  diferenca de tamanho: {len(after)-len(before):+d} chars")
    return 1


if __name__ == "__main__":
    sys.exit(main())
