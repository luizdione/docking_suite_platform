#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gerar_manual_pdf.py — Gera o manual_skills.pdf a partir de conteudo embutido.

Usa fpdf2 (pip install "fpdf2==2.5.7"). Tenta a fonte Unicode DejaVu (acentos,
setas, travessao); se nao achar, cai para Helvetica com sanitizacao latin-1.
Saida: manual_skills.pdf ao lado deste script (pasta .claude do repo).

Uso:
    pip install "fpdf2==2.5.7"
    python gerar_manual_pdf.py
"""
import os
import shutil
import tempfile

from fpdf import FPDF

AQUI = os.path.dirname(os.path.abspath(__file__))
SAIDA = os.path.join(AQUI, "manual_skills.pdf")

# --- Paleta -----------------------------------------------------------------
TEAL_DARK = (11, 85, 99)
TEAL = (11, 114, 133)
ALT = (242, 246, 247)
CODE_BG = (15, 43, 51)
CODE_FG = (232, 241, 242)
CALL_BG = (255, 248, 230)
CALL_BD = (224, 168, 0)
TEXT = (31, 41, 51)
MUTED = (102, 115, 128)
WHITE = (255, 255, 255)
LINE = (183, 207, 214)

# --- Fontes -----------------------------------------------------------------
DEJAVU_DIRS = [
    "/usr/share/fonts/truetype/dejavu",
    "C:/Windows/Fonts",
    "/Library/Fonts",
]
_FILES = {
    "": "DejaVuSans.ttf",
    "B": "DejaVuSans-Bold.ttf",
    "M": "DejaVuSansMono.ttf",
    "MB": "DejaVuSansMono-Bold.ttf",
}


def achar_fontes():
    for d in DEJAVU_DIRS:
        if all(os.path.exists(os.path.join(d, f)) for f in _FILES.values()):
            return d
    return None


class PDF(FPDF):
    def __init__(self, uni):
        super().__init__(orientation="P", unit="mm", format="A4")
        self.uni = uni
        self.set_auto_page_break(True, margin=16)
        self.set_margins(18, 16, 18)

    def footer(self):
        self.set_y(-12)
        self.set_font(self.fam, "", 8)
        self.set_text_color(*MUTED)
        self.cell(0, 6, self.S(f"Manual das Minhas Skills — pagina {self.page_no()}"),
                  0, 0, "C")

    # fam/famM definidos apos criar o objeto
    def S(self, txt):
        """Sanitiza para latin-1 quando estiver sem fonte Unicode."""
        if self.uni:
            return txt
        rep = {"—": "-", "–": "-", "→": "->", "←": "<-", "»": ">>", "“": '"',
               "”": '"', "‘": "'", "’": "'", "…": "...", "·": "-", "•": "-"}
        for k, v in rep.items():
            txt = txt.replace(k, v)
        return txt.encode("latin-1", "replace").decode("latin-1")


def build():
    fdir = achar_fontes()
    uni = fdir is not None
    pdf = PDF(uni)

    if uni:
        tmp = tempfile.mkdtemp()
        for st, fn in _FILES.items():
            shutil.copy(os.path.join(fdir, fn), os.path.join(tmp, fn))
        pdf.add_font("DV", "", os.path.join(tmp, _FILES[""]), uni=True)
        pdf.add_font("DV", "B", os.path.join(tmp, _FILES["B"]), uni=True)
        pdf.add_font("DVM", "", os.path.join(tmp, _FILES["M"]), uni=True)
        pdf.add_font("DVM", "B", os.path.join(tmp, _FILES["MB"]), uni=True)
        pdf.fam, pdf.famM = "DV", "DVM"
    else:
        pdf.fam, pdf.famM = "helvetica", "courier"

    pdf.alias_nb_pages()
    pdf.add_page()

    EPW = pdf.w - 2 * pdf.l_margin  # largura util

    # ---- helpers de layout ----
    def h1_bar(title, sub, sub2):
        pdf.set_fill_color(*TEAL_DARK)
        y0 = pdf.get_y()
        pdf.rect(pdf.l_margin, y0, EPW, 20, "F")
        pdf.set_xy(pdf.l_margin + 3, y0 + 3)
        pdf.set_font(pdf.fam, "B", 20)
        pdf.set_text_color(*WHITE)
        pdf.cell(0, 9, pdf.S(title), 0, 2)
        pdf.set_font(pdf.fam, "", 10)
        pdf.cell(0, 5, pdf.S(sub), 0, 2)
        pdf.set_y(y0 + 22)
        pdf.set_font(pdf.fam, "", 9)
        pdf.set_text_color(*MUTED)
        pdf.cell(0, 5, pdf.S(sub2), 0, 1)
        pdf.ln(3)

    def h2(txt):
        if pdf.get_y() > pdf.h - 40:
            pdf.add_page()
        pdf.ln(3)
        pdf.set_font(pdf.fam, "B", 14)
        pdf.set_text_color(*TEAL_DARK)
        pdf.cell(0, 8, pdf.S(txt), 0, 1)
        pdf.set_draw_color(*LINE)
        pdf.set_line_width(0.3)
        y = pdf.get_y()
        pdf.line(pdf.l_margin, y, pdf.l_margin + EPW, y)
        pdf.ln(2)

    def h3(txt):
        pdf.ln(1)
        pdf.set_font(pdf.fam, "B", 11.5)
        pdf.set_text_color(*TEAL)
        pdf.cell(0, 6, pdf.S(txt), 0, 1)

    def para(txt):
        pdf.set_font(pdf.fam, "", 10.5)
        pdf.set_text_color(*TEXT)
        pdf.multi_cell(EPW, 5.4, pdf.S(txt))
        pdf.ln(1)

    def bullets(items):
        pdf.set_font(pdf.fam, "", 10.5)
        pdf.set_text_color(*TEXT)
        for it in items:
            x = pdf.get_x()
            pdf.cell(5, 5.2, pdf.S("•"), 0, 0)
            pdf.multi_cell(EPW - 5, 5.2, pdf.S(it))
        pdf.ln(1)

    def code(lines):
        pdf.ln(1)
        n = len(lines)
        lh = 4.6
        pad = 2.2
        hbox = n * lh + 2 * pad
        if pdf.get_y() + hbox > pdf.h - pdf.b_margin:
            pdf.add_page()
        x, y = pdf.l_margin, pdf.get_y()
        pdf.set_fill_color(*CODE_BG)
        pdf.rect(x, y, EPW, hbox, "F")
        pdf.set_font(pdf.famM, "", 9)
        pdf.set_text_color(*CODE_FG)
        pdf.set_xy(x + 3, y + pad)
        for ln_ in lines:
            pdf.set_x(x + 3)
            pdf.cell(EPW - 6, lh, pdf.S(ln_), 0, 2)
        pdf.set_y(y + hbox)
        pdf.ln(2)

    def callout(txt):
        pdf.set_font(pdf.fam, "", 10)
        lines = pdf.multi_cell(EPW - 8, 5, pdf.S(txt), split_only=True)
        hbox = len(lines) * 5 + 5
        if pdf.get_y() + hbox > pdf.h - pdf.b_margin:
            pdf.add_page()
        x, y = pdf.l_margin, pdf.get_y()
        pdf.set_fill_color(*CALL_BG)
        pdf.rect(x, y, EPW, hbox, "F")
        pdf.set_fill_color(*CALL_BD)
        pdf.rect(x, y, 1.5, hbox, "F")
        pdf.set_xy(x + 5, y + 2.5)
        pdf.set_text_color(*TEXT)
        pdf.multi_cell(EPW - 8, 5, pdf.S(txt))
        pdf.set_y(y + hbox)
        pdf.ln(2)

    def table(headers, rows, widths):
        lh = 4.8
        pad = 1.4

        def draw_header():
            pdf.set_font(pdf.fam, "B", 9.5)
            pdf.set_fill_color(*TEAL)
            pdf.set_text_color(*WHITE)
            pdf.set_draw_color(*LINE)
            pdf.set_line_width(0.2)
            y = pdf.get_y()
            x = pdf.l_margin
            hh = lh + 2 * pad
            for w, htxt in zip(widths, headers):
                pdf.rect(x, y, w, hh, "DF")
                pdf.set_xy(x + 1.5, y + pad)
                pdf.multi_cell(w - 3, lh, pdf.S(htxt), 0, "L")
                x += w
                pdf.set_y(y)
            pdf.set_y(y + hh)

        draw_header()
        pdf.set_font(pdf.fam, "", 9.3)
        alt = False
        for row in rows:
            # calcula altura da linha
            cell_lines = []
            for w, cell in zip(widths, row):
                ls = pdf.multi_cell(w - 3, lh, pdf.S(str(cell)), split_only=True)
                cell_lines.append(ls)
            nlines = max(len(ls) for ls in cell_lines)
            rh = nlines * lh + 2 * pad
            if pdf.get_y() + rh > pdf.h - pdf.b_margin:
                pdf.add_page()
                draw_header()
                pdf.set_font(pdf.fam, "", 9.3)
            y = pdf.get_y()
            x = pdf.l_margin
            pdf.set_draw_color(*LINE)
            pdf.set_text_color(*TEXT)
            for w, ls in zip(widths, cell_lines):
                if alt:
                    pdf.set_fill_color(*ALT)
                    pdf.rect(x, y, w, rh, "DF")
                else:
                    pdf.set_fill_color(*WHITE)
                    pdf.rect(x, y, w, rh, "DF")
                pdf.set_xy(x + 1.5, y + pad)
                pdf.multi_cell(w - 3, lh, "\n".join(ls), 0, "L")
                x += w
                pdf.set_y(y)
            pdf.set_y(y + rh)
            alt = not alt
        pdf.ln(2)

    # ============ CONTEUDO ============
    h1_bar("Manual das Minhas Skills Pessoais",
           "Claude Code — Skills de memoria de sessao — NQ-Docking Suite",
           "Autor: Luiz Dione   |   Skills em: C:\\Users\\luizd\\.claude\\skills")

    h2("1. O que sao 'skills pessoais' e onde ficam")
    para("Uma skill e uma pastinha com um arquivo SKILL.md que ensina o Claude Code a "
         "executar uma tarefa. Quando a skill fica em C:\\Users\\luizd\\.claude\\skills\\ "
         "ela e pessoal (global): vale em TODOS os seus projetos, nao so na NQ-Docking Suite.")
    table(["Escopo", "Local", "Vale para"],
          [["Pessoal (global)", "C:\\Users\\luizd\\.claude\\skills\\<nome>\\SKILL.md",
            "Todos os seus projetos"],
           ["Projeto", "<repo>\\.claude\\skills\\<nome>\\SKILL.md",
            "So aquele projeto (versionado no git)"]],
          [34, 100, 40])
    para("Como chamar: digite '/' seguido do nome da skill. Ex.: /salvar_contexto. "
         "Voce tambem pode passar um foco: /salvar_contexto modulo ADMET.")

    h2("2. As duas skills e o que cada uma faz")
    h3("/salvar_contexto")
    para("Gera/atualiza a memoria da sessao atual em arquivos Markdown do projeto (objetivo, "
         "progresso, decisoes, arquivos tocados e pendencias). Serve para voce reduzir o "
         "contexto depois — com /compact, /clear ou /branch — sem perder a continuidade.")
    para("Quando usar: contexto grande (perto do limite / 1M no Opus), antes de "
         "compactar/limpar/forkar, ou ao encerrar um bloco de trabalho.")
    h3("/salvar_contexto_compactar")
    para("Faz o mesmo que a /salvar_contexto (grava a memoria) e, no fim, entrega o comando "
         "/compact pronto (com foco) para voce digitar — o fluxo 'salvar e compactar' de uma vez.")
    callout("Importante: nenhuma skill executa /compact sozinha. /compact, /clear e /branch "
            "sao comandos que VOCE digita no Claude Code. A skill prepara tudo (salva a "
            "memoria) e te diz o comando exato para digitar em seguida.")

    h2("3. Quais arquivos .md sao atualizados")
    para("As duas skills mexem nos mesmos tres arquivos, sempre dentro do projeto atual:")
    table(["Arquivo", "Papel", "Como e atualizado"],
          [["docs/historico/CONTEXTO_SESSAO.md",
            "Arquivo vivo — estado mais recente (o que se rele ao retomar)",
            "Sobrescrito a cada execucao"],
           ["docs/historico/BRIEFING_<data>_<tag>.md",
            "Snapshot datado — historico imutavel da sessao",
            "Novo a cada vez; nunca sobrescreve os anteriores"],
           ["claude.md (raiz)",
            "Ponteiro de compactacao",
            "Ganha a secao '## Instrucoes de Compactacao' so na 1a vez (idempotente)"]],
          [56, 66, 52])
    para("» Para retomar do zero, o arquivo-chave e docs/historico/CONTEXTO_SESSAO.md.")

    h2("4. O hook de auto-load (memoria automatica ao abrir a sessao)")
    para("O projeto tem um hook SessionStart (.claude/settings.json + "
         ".claude/hooks/carregar_contexto.py) que injeta automaticamente o conteudo de "
         "docs/historico/CONTEXTO_SESSAO.md no inicio de TODA sessao nova, retomada, apos "
         "/clear e apos /compact. Se o arquivo nao existir, o hook nao faz nada (silencioso). "
         "Assim, ao abrir a sessao, o Claude ja 'lembra' onde voce parou — sem digitar nada.")
    para("Teto de seguranca: o hook injeta no maximo ~20 mil caracteres, para nao inflar o "
         "contexto logo na abertura. Mantenha o CONTEXTO_SESSAO.md enxuto.")

    h2("5. Fluxo completo recomendado")
    table(["#", "Situacao", "O que fazer"],
          [["1", "Contexto ficando grande (perto de 1M)",
            "/salvar_contexto_compactar (ou /salvar_contexto)"],
           ["2", "Reduzir de fato o contexto",
            "Digitar /compact foque em <tema> (a skill da o comando pronto)"],
           ["3", "(Opcional) fazer a memoria viajar",
            "git add/commit dos .md em docs/historico/"],
           ["4", "Abrir sessao nova do zero",
            "O hook carrega a memoria sozinho. Detalhe: 'releia CONTEXTO_SESSAO.md'"]],
          [10, 70, 94])

    h2("6. Instalar / atualizar as skills pessoais")
    para("As skills sao desenvolvidas e versionadas no repositorio e depois copiadas para a "
         "sua pasta pessoal. Rode o instalador (no Windows) a partir da pasta .claude do repo:")
    code(['cd "C:\\Users\\luizd\\OneDrive\\Data_Science\\Projetos\\docking\\.claude"',
          '.\\instalar_skills_pessoais.ps1'])
    para("O instalador copia salvar_contexto e salvar_contexto_compactar para "
         "C:\\Users\\luizd\\.claude\\skills\\ e copia este manual para "
         "...\\Projetos\\manual_skills.pdf. Para reinstalar depois de editar uma skill, "
         "e so rodar de novo (ele sobrescreve).")
    callout("Como estao em escopo pessoal, as skills passam a valer em qualquer projeto seu. "
            "Dentro da NQ-Docking Suite elas tambem existem no escopo do projeto (repo) — "
            "funciona igual; se preferir, pode manter so a versao pessoal.")

    h2("7. Atalho rapido para a pasta da skill")
    para("Ha um script que cria um atalho do Windows apontando para a pasta da skill, dentro "
         "de Projetos:")
    code(['cd "...\\docking\\.claude\\skills\\salvar_contexto"',
          '.\\criar_atalho.ps1   # cria "salvar_contexto (skill).lnk" em ...\\Projetos'])

    h2("8. Referencia rapida (cola)")
    table(["Comando", "Efeito"],
          [["/salvar_contexto", "Grava a memoria nos 3 arquivos .md"],
           ["/salvar_contexto_compactar", "Grava a memoria e entrega o /compact pronto"],
           ["/compact foque em X", "Reduz o contexto mantendo o essencial (voce digita)"],
           ["/clear  ->  /resume", "Zera o contexto; depois retoma a sessao anterior"],
           ["/branch", "Forka a sessao (copia com novo ID), mantem a original"],
           ["/context", "Mostra o que esta ocupando o contexto"]],
          [60, 114])

    pdf.output(SAIDA)
    print("PDF gerado:", SAIDA, "| unicode:", uni)


if __name__ == "__main__":
    build()
