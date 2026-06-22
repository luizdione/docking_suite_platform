#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pipeline_svg.py — gera o diagrama SVG do pipeline de cada subprojeto da plataforma
DOCKING (padrao guarda-chuva), espelhando o gerador irmao da plataforma epigenoma
(epigenoma/scripts/python/pipeline_svg.py).

Le o registry `subprojeto/subprojeto_config.json` e o `subprojeto/<ID>/status.json`
de cada subprojeto e pinta cada etapa por status (done/pending/blocked/data).

Saidas (por subprojeto ID, ex. SUB01):
    docking_project/figures/pipeline_<id_minusculo>.svg   (servido pelo dashboard)
    subprojeto/<ID>/pipeline_<id_minusculo>.svg           (copia p/ submissao)

Sem dependencias externas: escreve SVG como texto. Reexecutavel a cada versionamento.
Uso:
    python subprojeto/pipeline_svg.py            # todos os subprojetos
    python subprojeto/pipeline_svg.py SUB01      # apenas um
"""
import json
import sys
from pathlib import Path
from html import escape

BASE = Path(__file__).resolve().parents[1]              # .../docking
CONF = BASE / "subprojeto" / "subprojeto_config.json"

ST = {
    "done":    {"bg": "#d4edda", "br": "#1e8449", "tag": "concluido"},
    "pending": {"bg": "#fff3cd", "br": "#b9770e", "tag": "a fazer"},
    "blocked": {"bg": "#e2e3e5", "br": "#6c757d", "tag": "bloqueado / futuro"},
    "data":    {"bg": "#cfe2f3", "br": "#1f5fa8", "tag": "dados"},
}
W, H = 1500, 760


def anchor(box, side):
    x, y, w, h = box["x"], box["y"], box["w"], box["h"]
    return {"n": (x + w / 2, y), "s": (x + w / 2, y + h),
            "e": (x + w, y + h / 2), "w": (x, y + h / 2)}[side]


def svg_box(b):
    c = ST.get(b["status"], ST["pending"])
    x, y, w, h = b["x"], b["y"], b["w"], b["h"]
    return (
        f'<g><rect x="{x}" y="{y}" rx="9" ry="9" width="{w}" height="{h}" '
        f'fill="{c["bg"]}" stroke="{c["br"]}" stroke-width="2"/>'
        f'<text x="{x + w/2}" y="{y + 22}" text-anchor="middle" '
        f'font-family="Segoe UI, Arial" font-size="14" font-weight="700" '
        f'fill="#1a1a1a">{escape(b["title"])}</text>'
        f'<text x="{x + w/2}" y="{y + 40}" text-anchor="middle" '
        f'font-family="Segoe UI, Arial" font-size="11" fill="#444">{escape(b["sub"])}</text></g>'
    )


def svg_arrow(boxes, a, sa, b, sb):
    x1, y1 = anchor(boxes[a], sa)
    x2, y2 = anchor(boxes[b], sb)
    return (f'<path d="M {x1} {y1} L {x2} {y2}" fill="none" stroke="#7a7a7a" '
            f'stroke-width="2" marker-end="url(#ah)"/>')


def lane_label(y, text, color):
    return (f'<text x="20" y="{y}" font-family="Segoe UI, Arial" font-size="13" '
            f'font-weight="700" fill="{color}" transform="rotate(-90 20 {y})">'
            f'{escape(text)}</text>')


def status_of(steps, key, default="pending"):
    for s in steps:
        if s.get("etapa") == key:
            return s.get("status", default)
    return default


def build(sp_id, sp, steps):
    nome = sp.get("nome_especifico", sp_id)
    titulo = sp.get("titulo_paper_atual", nome)

    def st(k):
        return status_of(steps, k)

    # Lanes: DADOS / DOCKING / MD (n=3) / ANALISE / PUBLICAR
    boxes = {
        "d_lig":  {"x": 40,  "y": 80,  "w": 250, "h": 60, "title": "2-SCN-NQ", "sub": "C11H5NO2S - MW 215,23", "status": "data"},
        "d_rec":  {"x": 320, "y": 80,  "w": 300, "h": 60, "title": "Receptores (9 alvos)", "sub": "TcTR, LiTR, hGR, hCatL...", "status": "data"},
        "d_ff":   {"x": 650, "y": 80,  "w": 250, "h": 60, "title": "FF FAD/LIG", "sub": "smoke_test_fad XML", "status": "data"},

        "dock":   {"x": 40,  "y": 210, "w": 280, "h": 64, "title": "Docking Vina 1.2.7", "sub": "painel 9 alvos", "status": st("docking_painel_9_alvos")},

        "md_tc":  {"x": 40,  "y": 330, "w": 230, "h": 60, "title": "MD TcTR (n=3)", "sub": "50 ns - 2,26 A FAD", "status": st("md_subversivo_n3_TcTR")},
        "md_li":  {"x": 290, "y": 330, "w": 230, "h": 60, "title": "MD LiTR (n=3)", "sub": "50 ns - 2,43 A FAD", "status": st("md_subversivo_n3_LiTR")},
        "md_hg":  {"x": 540, "y": 330, "w": 230, "h": 60, "title": "MD hGR (n=3)", "sub": "50 ns - 9,28 A FAD", "status": st("md_subversivo_n3_hGR")},
        "md_inh": {"x": 790, "y": 330, "w": 230, "h": 60, "title": "MD inibidor (n=3)", "sub": "comparacao A x B", "status": st("md_inibidor_n3")},

        "agg":    {"x": 40,  "y": 450, "w": 270, "h": 60, "title": "Agregacao n=3", "sub": "n3_resultados.json", "status": st("agregacao_n3")},
        "sel":    {"x": 330, "y": 450, "w": 290, "h": 60, "title": "Seletividade lig-FAD", "sub": "parasito ~2,3 A vs humano ~9,3 A", "status": st("metrica_seletividade_lig_FAD")},
        "qm":     {"x": 640, "y": 450, "w": 270, "h": 60, "title": "QM/MM ou dG", "sub": "estado radicalar / energia", "status": st("qmmm_ou_dG_ligacao")},

        "paper":  {"x": 40,  "y": 580, "w": 300, "h": 60, "title": "Manuscrito (build_paper)", "sub": "ACS Infect. Dis. - EN/PT", "status": st("manuscrito_build_paper")},
        "figs":   {"x": 360, "y": 580, "w": 270, "h": 60, "title": "Figuras 300 dpi + ESM", "sub": "submissao", "status": st("figuras_300dpi_e_suplementar")},
        "cover":  {"x": 650, "y": 580, "w": 270, "h": 60, "title": "Cover + declaracoes", "sub": "ORCID/COI/CRediT", "status": st("cover_letter_declaracoes")},
    }
    arrows = [
        ("d_lig", "s", "dock", "n"), ("d_rec", "s", "dock", "n"), ("d_ff", "s", "md_tc", "n"),
        ("dock", "s", "md_tc", "n"), ("dock", "s", "md_li", "n"), ("dock", "s", "md_hg", "n"),
        ("md_tc", "s", "agg", "n"), ("md_li", "s", "agg", "n"), ("md_hg", "s", "agg", "n"),
        ("agg", "e", "sel", "w"), ("sel", "e", "qm", "w"),
        ("sel", "s", "paper", "n"), ("paper", "e", "figs", "w"), ("figs", "e", "cover", "w"),
    ]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Segoe UI, Arial">',
        '<defs><marker id="ah" markerWidth="10" markerHeight="10" refX="8" refY="3" '
        'orient="auto" markerUnits="strokeWidth"><path d="M0,0 L8,3 L0,6 z" fill="#7a7a7a"/></marker></defs>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>',
        f'<text x="{W/2}" y="34" text-anchor="middle" font-size="20" font-weight="800" fill="#15406b">'
        f'{escape(sp_id)} - pipeline e estado para publicacao</text>',
        f'<text x="{W/2}" y="56" text-anchor="middle" font-size="12" fill="#666">'
        f'{escape(titulo[:120])}</text>',
        lane_label(110, "DADOS", "#1f5fa8"),
        lane_label(242, "DOCKING", "#1e8449"),
        lane_label(360, "MD (n=3)", "#1e8449"),
        lane_label(480, "ANALISE", "#1e8449"),
        lane_label(610, "PUBLICAR", "#b9770e"),
    ]
    for a, sa, b, sb in arrows:
        parts.append(svg_arrow(boxes, a, sa, b, sb))
    for bid in boxes:
        parts.append(svg_box(boxes[bid]))
    # legenda
    lx, ly = 40, 710
    parts.append(f'<text x="{lx}" y="{ly-8}" font-size="13" font-weight="700" fill="#333">Legenda:</text>')
    for i, k in enumerate(["done", "pending", "blocked", "data"]):
        c = ST[k]; bx = lx + 90 + i * 300
        parts.append(f'<rect x="{bx}" y="{ly-20}" width="20" height="14" rx="3" fill="{c["bg"]}" '
                     f'stroke="{c["br"]}" stroke-width="2"/>')
        parts.append(f'<text x="{bx+28}" y="{ly-9}" font-size="12" fill="#333">{escape(c["tag"])}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def main():
    conf = json.loads(CONF.read_text(encoding="utf-8"))
    only = sys.argv[1] if len(sys.argv) > 1 else None
    figs = BASE / "docking_project" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    for sp_id, sp in conf.get("subprojetos", {}).items():
        if only and sp_id != only:
            continue
        st_path = BASE / "subprojeto" / sp_id / "status.json"
        steps = []
        if st_path.exists():
            steps = json.loads(st_path.read_text(encoding="utf-8")).get("etapas", [])
        svg = build(sp_id, sp, steps)
        out1 = figs / f"pipeline_{sp_id.lower()}.svg"
        out2 = BASE / "subprojeto" / sp_id / f"pipeline_{sp_id.lower()}.svg"
        out1.write_text(svg, encoding="utf-8")
        out2.write_text(svg, encoding="utf-8")
        print("SVG salvo:", out1)
        print("SVG salvo:", out2)


if __name__ == "__main__":
    main()
