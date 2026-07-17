#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""export_web.py - empacota o BUNDLE ESTATICO publicavel (NQ-Docking Viewer).

Gera um site 100% estatico (HTML/JS/JSON) publicavel na Hostinger
(`docking.ifrj-crj-geneticamolecular.online` ou subpasta `/docking/`), em modo
**somente leitura** + calculadora ADMET client-side, sem Python/GPU/rede no servidor.

FIDELIDADE AO PAPER (2026-07-03): os dados cientificos (docking, MD, ADMET, seletividade
e as FIGURAS) sao os do **manuscrito ACS submetido** (`2scnnq_manuscript_EN.tex` /
`2scnnq_acs_submission.pdf`), NAO os dados de trabalho do dashboard. Os valores abaixo
foram extraidos e conferidos contra o .tex submetido (Tabela tab:dock e secoes
Docking/MD/ADMET). Metadados nao-cientificos (biblioteca de compostos, subprojetos, docs
da Ajuda) continuam vindo do dashboard via test_client.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
VIEWER_DIR = SCRIPT_DIR / "web" / "viewer"

# --- Projeto do paper (fonte de verdade cientifica) --------------------------
PAPER_DIR = Path(os.environ.get("NQ_PAPER_DIR", r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-SNC-NQ"))
PAPER_FIGS_DIR = PAPER_DIR / "PAPER" / "figures"
SUBMISSION_PDF_CANDIDATES = [
    Path(r"C:\Users\luizd\OneDrive\Data_Science\Projetos\artigos_submissao\naftoquinona_2SCN-NQ\2scnnq_acs_submission.pdf"),
    PAPER_DIR / "PAPER" / "2scnnq_acs_submission.pdf",
]

# =============================================================================
# DADOS CANONICOS DO MANUSCRITO SUBMETIDO (conferidos contra 2scnnq_manuscript_EN.tex)
# =============================================================================
PAPER = {
    "compound": {
        "name": "2-SCN-NQ",
        "full": "2-Thiocyanato-1,4-naphthoquinone",
        "formula": "C11H5NO2S",
        "mw": 215.23,
        "smiles": "O=C1C=C(SC#N)C(=O)c2ccccc21",
    },
    # Tabela tab:dock (manuscrito submetido), ordenada por dG. Ki em uM.
    # rmsd2 = Top-2 pose RMSD (A); None quando so 1 pose (hDHODH).
    "docking": [
        {"code": "LiTR", "label": "TR (LiTR)", "organism": "L. infantum", "pdb": "2JK6", "dg": -7.55, "ki_uM": 2.90, "rmsd2": 22.9, "category": "leishmania", "cofactors": ["FAD"]},
        {"code": "hGR", "label": "GR (hGR)", "organism": "H. sapiens", "pdb": "3DJJ", "dg": -7.16, "ki_uM": 5.64, "rmsd2": 4.8, "category": "human", "cofactors": ["FAD"]},
        {"code": "TcTR", "label": "TR (TcTR)", "organism": "T. cruzi", "pdb": "1BZL", "dg": -6.86, "ki_uM": 9.42, "rmsd2": 2.7, "category": "parasite", "cofactors": ["FAD"]},
        {"code": "hCatL", "label": "Cathepsin L", "organism": "H. sapiens", "pdb": "2XU3", "dg": -6.68, "ki_uM": 12.74, "rmsd2": 14.0, "category": "human", "cofactors": []},
        {"code": "TcFeSOD", "label": "Fe-SOD (TcFeSOD)", "organism": "T. cruzi", "pdb": "4DVH", "dg": -6.26, "ki_uM": 25.89, "rmsd2": 4.6, "category": "parasite", "cofactors": ["FE"]},
        {"code": "hDHODH", "label": "DHODH (hDHODH)", "organism": "H. sapiens", "pdb": "1D3H", "dg": -6.02, "ki_uM": 38.89, "rmsd2": None, "category": "human", "cofactors": ["FMN"]},
        {"code": "TcCruzain", "label": "Cruzain", "organism": "T. cruzi", "pdb": "3KKU", "dg": -6.00, "ki_uM": 39.82, "rmsd2": 13.6, "category": "parasite", "cofactors": []},
        {"code": "LdDHODH", "label": "DHODH (LdDHODH)", "organism": "L. donovani", "pdb": "3C61", "dg": -5.29, "ki_uM": 130.9, "rmsd2": 4.7, "category": "leishmania", "cofactors": ["FMN"]},
        {"code": "hMnSOD", "label": "Mn-SOD (hMnSOD)", "organism": "H. sapiens", "pdb": "1LUV", "dg": -4.85, "ki_uM": 280.2, "rmsd2": 2.8, "category": "human", "cofactors": ["MN"]},
    ],
    # Dinamica molecular n=3 (secoes MD do manuscrito). Distancias em A.
    "md": {
        "subversive": [
            {"target": "TcTR", "fad": 2.26, "sd": 0.07, "dithiol_pct": 97, "category": "parasite"},
            {"target": "LiTR", "fad": 2.43, "sd": 0.23, "dithiol_pct": 95, "category": "leishmania"},
            {"target": "hGR", "fad": 9.28, "sd": 1.44, "dithiol_pct": 0, "category": "human"},
        ],
        "inhibitor": [
            {"target": "TcTR", "fad": 13.7, "dithiol": 16.7},
            {"target": "LiTR", "fad": 11.2},
        ],
        "stats": "Welch ANOVA p = 0,012; hGR vs. cada TR parasitaria p = 0,038 (Welch t, Holm); as duas TRs entre si p = 0,33; Cohen d ~ 7.",
        "criterion": "Redox-competente: media ligante-FAD < 3 A E contato com o ditiol catalitico em > 50% dos frames. So as poses subversivas das TRs parasitarias qualificam.",
        "protocol": "AMBER ff14SB + TIP3P, 150 mM NaCl, 300 K; OpenMM 8.3 (HMR, passo 4 fs); 50 ns por corrida; subversivo em triplicata (n=3), inibidor n=1; 12 trajetorias.",
        "contacts_tctr": [["Cys52", 2.5], ["Cys57", 2.5], ["Val58", None], ["Ser14", None], ["Gly56", None], ["Thr51", None], ["Pro59", None], ["Val55", None]],
        "fold_note": "RMSD do backbone por cadeia ~1,4 A (fold preservado); RMSD intramolecular do ligante < 0,6 A (pose estavel).",
    },
    # ADMET in silico (secao ADMET do manuscrito).
    "admet": {
        "MW": 215.23, "cLogP": 2.16, "TPSA": 57.9, "HBD": 0, "HBA": 4, "RB": 1, "QED": 0.935,
        "lipinski_viol": 0, "veber": True, "egan": True, "muegge": True, "ghose": False,
        "bioavailability": 1.00, "hia": 54.3, "logS": -3.21, "logBB": -0.89, "ppb": 75.4, "vd": 10.1,
        "herg": "baixa", "ames": False, "cyp_risk": "alto", "pains": 3, "t_half": 48,
        "caco2": "alta", "bbb": True,
        "cyp_motifs": ["naftoquinona", "tiocianato", "tioeter"],
    },
    "selectivity_cellular_max": 36.6,
    "ic50_note": "Concordancia com a potencia celular: IC50 5,64 uM (L. amazonensis) e 10,6 uM (T. cruzi).",
}

# Figuras OFICIAIS do paper (arquivo -> rotulo). Atencao: no manuscrito a numeracao
# de arquivo nao segue a ordem impressa (fig11=docking, fig12=ADMET, fig10=MD).
PAPER_FIGS = {
    "structure": ("2scnnq_structure.svg", "Estrutura 2D de 2-SCN-NQ"),
    "docking": ("fig11.png", "Docking: ranking de afinidade (dG/Ki), indice de seletividade e convergencia de poses"),
    "admet": ("fig12.png", "Perfil ADMET in silico (filtros de drug-likeness, absorcao/distribuicao, seguranca)"),
    "md_tctr": ("fig10.png", "MD do complexo TcTR-2-SCN-NQ: RMSD e metricas do modo de ligacao"),
    "md_matrix": ("fig_md_matrix.png", "Seletividade dinamica parasita-vs-humano: matriz de competencia redox e distancia ligante-FAD por sistema"),
    "interaction_2d": ("fig13.png", "Mapas de interacao 3D no sitio da TcTR: modo inibidor (periferico) vs. subversivo (face si do FAD)"),
    "interaction_3d": ("fig_md3d.png", "Mapas de interacao 3D dos seis sistemas de MD (enzima x modo)"),
    "graphical_abstract": ("graphical_abstract.png", "Graphical abstract"),
}

# -----------------------------------------------------------------------------
DROP_KEYS = {"_dir", "_has_mol2", "_has_pdbqt", "path", "docking_dir", "figures_dir",
             "logs_dir", "results_file", "texts_dir", "kit_dir"}
_ABS_PATH = re.compile(r"^[A-Za-z]:[\\/]|^/(home|Users|mnt|data)/")


def _looks_like_local_path(s):
    return bool(_ABS_PATH.match(s)) or ("\\" in s and ":" in s)


def sanitize(obj):
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items() if k not in DROP_KEYS}
    if isinstance(obj, list):
        return [sanitize(x) for x in obj]
    if isinstance(obj, str) and _looks_like_local_path(obj):
        return None
    return obj


def _write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def _get_json(client, url):
    r = client.get(url)
    if r.status_code != 200:
        return None
    try:
        return json.loads(r.get_data(as_text=True))
    except Exception:
        return None


def _copy_paper_figures(out):
    """Copia as figuras OFICIAIS do paper para figures/ e devolve o mapa rotulo->rel/legenda."""
    figs_out = out / "figures"
    figs_out.mkdir(parents=True, exist_ok=True)
    result = {}
    for key, (fname, caption) in PAPER_FIGS.items():
        src = PAPER_FIGS_DIR / fname
        if src.exists():
            shutil.copyfile(src, figs_out / fname)
            result[key] = {"rel": f"figures/{fname}", "caption": caption}
    return result


def export_bundle(out_dir, compound="2-SCN-NQ"):
    out = Path(out_dir)
    if not VIEWER_DIR.exists():
        raise FileNotFoundError(f"viewer nao encontrado em {VIEWER_DIR}")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(VIEWER_DIR, out, dirs_exist_ok=True)

    summary = {"out_dir": str(out), "figures": 0, "docs": 0, "pdf": False}

    # 1) DADOS CIENTIFICOS = manuscrito submetido -----------------------------
    _write_json(out / "data" / "docking.json", {
        "targets": PAPER["docking"],
        "si_note": "Indice de seletividade in silico = Ki(humano) / Ki(parasita). "
                   "A seletividade NAO decorre do docking energetico (hGR pontua entre as duas TRs); "
                   "e uma propriedade estrutural do canal de substrato da TR, confirmada pela MD.",
    })
    _write_json(out / "data" / "md.json", PAPER["md"])
    _write_json(out / "data" / "admet.json", {"admet": PAPER["admet"]})

    # 2) figuras oficiais do paper --------------------------------------------
    figs_map = _copy_paper_figures(out)
    summary["figures"] = len(figs_map)

    # 3) PDF de submissao: NAO incluir no viewer (o manuscrito nao e distribuido pela plataforma)
    pdf_rel = None

    # 4) metadados nao-cientificos do dashboard (compostos, subprojetos, docs) --
    sys.path.insert(0, str(SCRIPT_DIR))
    import dashboard_v2 as D
    client = D.app.test_client()

    comp = _get_json(client, "/compounds") or {}
    _write_json(out / "data" / "compounds.json", sanitize(comp))

    sp = _get_json(client, "/api/subprojects") or {}
    _write_json(out / "data" / "subprojects.json", sanitize(sp))
    active_sid = (sp.get("subprojetos") or [{}])[0].get("id")
    if active_sid:
        r = client.get(f"/api/subproject/{active_sid}/pipeline.svg")
        if r.status_code == 200:
            (out / "assets").mkdir(parents=True, exist_ok=True)
            (out / "assets" / "pipeline_SUB01.svg").write_bytes(r.get_data())

    docs_index = []
    (out / "docs").mkdir(parents=True, exist_ok=True)
    for d in (_get_json(client, "/help_list") or {}).get("docs", []):
        fname = d.get("file")
        if not fname:
            continue
        key = fname[:-3] if fname.lower().endswith(".md") else fname
        doc = _get_json(client, f"/help_doc/{key}")
        raw = (doc or {}).get("raw", "")
        if not raw:
            continue
        (out / "docs" / fname).write_text(raw, encoding="utf-8")
        docs_index.append({"id": key, "title": d.get("title", key), "file": fname})
    _write_json(out / "docs" / "index.json", {"docs": docs_index})
    summary["docs"] = len(docs_index)

    # 5) config.json -----------------------------------------------------------
    version = "0.0.0"
    vfile = SCRIPT_DIR / "VERSION"
    if vfile.exists():
        version = vfile.read_text(encoding="utf-8").strip()
    subproject = (sp.get("subprojetos") or [{}])[0]
    config = {
        "title": "Docking Suite IFRJ · 2-SCN-NQ, substrato subversivo da tripanotiona redutase",
        "version": version,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "source": "Resultados in silico (docking, dinâmica molecular e ADMET) do 2-SCN-NQ.",
        "compound": PAPER["compound"],
        "paper_title": subproject.get("titulo_paper"),
        "periodico_alvo": subproject.get("periodico_alvo"),
        "selectivity_cellular_max": PAPER["selectivity_cellular_max"],
        "ic50_note": PAPER["ic50_note"],
        "figs": figs_map,
        "submission_pdf": pdf_rel,
        "subproject": {
            "id": subproject.get("id"),
            "titulo_paper": subproject.get("titulo_paper"),
            "periodico_alvo": subproject.get("periodico_alvo"),
            "fase_atual": subproject.get("fase_atual"),
            "status": subproject.get("status"),
        },
    }
    _write_json(out / "config.json", config)
    return summary


def main():
    ap = argparse.ArgumentParser(description="Empacota o NQ-Docking Viewer estatico (fiel ao paper).")
    ap.add_argument("--out", default="dist_web")
    ap.add_argument("--compound", default="2-SCN-NQ")
    args = ap.parse_args()
    s = export_bundle(args.out, compound=args.compound)
    print("[OK] bundle gerado (fiel ao manuscrito submetido):")
    for k, v in s.items():
        print(f"     {k}: {v}")
    print(f"\nEnvie o conteudo de {s['out_dir']}/ para o Hostinger (ver README_DEPLOY.md).")


if __name__ == "__main__":
    main()
