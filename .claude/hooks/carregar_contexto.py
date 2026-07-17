#!/usr/bin/env python3
"""
carregar_contexto.py — Hook SessionStart do Claude Code.

Injeta a memoria da sessao anterior no contexto ao iniciar/retomar/limpar/compactar
a sessao. Le docs/historico/CONTEXTO_SESSAO.md (relativo a raiz do projeto = cwd do
hook) e imprime no stdout; o Claude Code adiciona a saida ao contexto inicial.

Silencioso (nao imprime nada, exit 0) se o arquivo nao existir ou estiver vazio.
Limita o tamanho injetado para nao inflar o contexto logo na abertura.
"""
import os
import sys

CAMINHO = os.path.join("docs", "historico", "CONTEXTO_SESSAO.md")
MAX_CHARS = 20000  # teto de seguranca para nao inflar o contexto inicial


def main() -> int:
    if not os.path.exists(CAMINHO):
        return 0
    try:
        with open(CAMINHO, encoding="utf-8") as fh:
            conteudo = fh.read().strip()
    except OSError:
        return 0
    if not conteudo:
        return 0

    if len(conteudo) > MAX_CHARS:
        conteudo = (
            conteudo[:MAX_CHARS]
            + "\n\n... [truncado — arquivo completo em docs/historico/CONTEXTO_SESSAO.md]"
        )

    print("## Memoria da sessao anterior (auto-carregada de docs/historico/CONTEXTO_SESSAO.md)")
    print()
    print(conteudo)
    print()
    print("> Continue pela secao 'Pendencias / proximos passos'. Para o detalhe completo,")
    print("> leia docs/historico/CONTEXTO_SESSAO.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
