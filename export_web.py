#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""export_web.py - empacota o BUNDLE ESTATICO publicavel (NQ-Docking Viewer).

Gera um site 100% estatico (HTML/JS/JSON) publicavel na Hostinger
(`docking.ifrj-crj-geneticamolecular.online` ou subpasta `/docking/`), em modo
**somente leitura** + calculadora ADMET client-side, sem Python/GPU/rede no servidor.

FIDELIDADE AO PAPER (2026-07-03): os dados cientificos (docking, MD, ADMET, seletividade
e as FIGURAS) sao os do **manuscrito ACS submetido** (`2scnnq_manuscript_EN.tex` /
`2scnnq_acs_submission.pdf`), NAO os dados de trabalho do dashboard.

T-0016 (Fase 8.2-8.3): esses dados deixaram de ser constantes hard-coded neste modulo
e passaram a viver CONGELADOS em `registry/campaigns/<id>.paper.json` (frozen: true),
separados dos resultados vivos. `export_bundle()` virou `build_viewer(campaign_id)`,
generica: carrega o paper congelado da campanha, ou exporta SEM secao de paper quando
a campanha ainda nao tem manuscrito (ex.: dados de trabalho). Metadados nao-cientificos
(biblioteca de compostos, subprojetos, docs da Ajuda) continuam vindo do dashboard via
test_client.
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
# Default historico do SUB01; build_viewer() resolve por campanha (ver _paper_figs_dir).
PAPER_DIR = Path(os.environ.get("NQ_PAPER_DIR", r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-SNC-NQ"))
PAPER_FIGS_DIR = PAPER_DIR / "PAPER" / "figures"
SUBMISSION_PDF_CANDIDATES = [
    Path(r"C:\Users\luizd\OneDrive\Data_Science\Projetos\artigos_submissao\naftoquinona_2SCN-NQ\2scnnq_acs_submission.pdf"),
    PAPER_DIR / "PAPER" / "2scnnq_acs_submission.pdf",
]


# =============================================================================
# DADOS CANONICOS DO MANUSCRITO -- agora CONGELADOS em JSON (T-0016 Fase 8.2)
# =============================================================================
def _load_paper(campaign_id):
    """Carrega registry/campaigns/<id>.paper.json (frozen), ou None se nao existir."""
    p = SCRIPT_DIR / "registry" / "campaigns" / f"{campaign_id}.paper.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


# Backward-compat: `export_web.PAPER`/`PAPER_FIGS` (usados por
# tools/verify_paper_vs_pipeline.py) agora sao servidos a partir do JSON congelado
# do SUB01, nao mais literais hard-coded aqui.
_SUB01_PAPER = _load_paper("SUB01") or {}
PAPER = {
    "compound": _SUB01_PAPER.get("compound"),
    "docking": _SUB01_PAPER.get("docking"),
    "md": _SUB01_PAPER.get("md"),
    "admet": _SUB01_PAPER.get("admet"),
    "selectivity_cellular_max": _SUB01_PAPER.get("selectivity_cellular_max"),
    "ic50_note": _SUB01_PAPER.get("ic50_note"),
}
# figuras oficiais: {chave: (arquivo, legenda)}
PAPER_FIGS = {k: tuple(v) for k, v in (_SUB01_PAPER.get("figs") or {}).items()}

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


def _paper_figs_dir(campaign_id, paper):
    """Diretorio das figuras oficiais do paper de uma campanha.

    Precedencia (preserva o comportamento historico do SUB01):
      env NQ_PAPER_DIR  ->  campaign.paths.paper_root  ->  default do modulo.
    """
    env = os.environ.get("NQ_PAPER_DIR")
    if env:
        return Path(env) / "PAPER" / "figures"
    try:
        import registry as _reg
        camp = _reg.get_campaign(campaign_id)
        proot = (camp.get("paths") or {}).get("paper_root")
        if proot:
            return Path(proot) / "figures"
    except Exception:
        pass
    return PAPER_FIGS_DIR


def _copy_paper_figures(out, campaign_id, paper):
    """Copia as figuras OFICIAIS do paper para figures/ e devolve o mapa rotulo->rel/legenda."""
    figs_out = out / "figures"
    figs_out.mkdir(parents=True, exist_ok=True)
    figs_dir = _paper_figs_dir(campaign_id, paper)
    result = {}
    for key, val in (paper.get("figs") or {}).items():
        fname, caption = val[0], val[1]
        src = figs_dir / fname
        if src.exists():
            shutil.copyfile(src, figs_out / fname)
            result[key] = {"rel": f"figures/{fname}", "caption": caption}
    return result


def build_viewer(campaign_id="SUB01", out_dir="dist_web", include_pipeline=True):
    """Empacota o viewer estatico de UMA campanha (T-0016 Fase 8.2-8.3).

    Se existir `registry/campaigns/<id>.paper.json` (frozen), usa os dados
    cientificos congelados do manuscrito. Caso contrario, exporta SEM secao de
    paper (dados de trabalho; abas de docking/MD/ADMET ausentes).
    """
    out = Path(out_dir)
    if not VIEWER_DIR.exists():
        raise FileNotFoundError(f"viewer nao encontrado em {VIEWER_DIR}")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(VIEWER_DIR, out, dirs_exist_ok=True)

    paper = _load_paper(campaign_id)
    summary = {"out_dir": str(out), "campaign": campaign_id,
               "has_paper": bool(paper), "figures": 0, "docs": 0, "pdf": False}

    # 1) DADOS CIENTIFICOS = manuscrito congelado (so quando ha paper) ---------
    figs_map = {}
    if paper:
        _write_json(out / "data" / "docking.json", {
            "targets": paper["docking"],
            "si_note": paper.get("docking_si_note", ""),
        })
        _write_json(out / "data" / "md.json", paper["md"])
        _write_json(out / "data" / "admet.json", {"admet": paper["admet"]})
        # 2) figuras oficiais do paper ----------------------------------------
        figs_map = _copy_paper_figures(out, campaign_id, paper)
        summary["figures"] = len(figs_map)

    # 3) PDF de submissao: NAO incluir no viewer (manuscrito nao e distribuido)
    pdf_rel = None

    # 4) metadados nao-cientificos do dashboard (compostos, subprojetos, docs) --
    sys.path.insert(0, str(SCRIPT_DIR))
    import dashboard_v2 as D
    client = D.app.test_client()

    comp = _get_json(client, f"/compounds?campaign={campaign_id}") or {}
    _write_json(out / "data" / "compounds.json", sanitize(comp))

    sp = _get_json(client, "/api/subprojects") or {}
    _write_json(out / "data" / "subprojects.json", sanitize(sp))
    # Subprojeto correspondente A campanha exportada (corrige o antigo bug do
    # indice fixo [0], que assumia o primeiro subprojeto da lista).
    subs = sp.get("subprojetos") or []
    # Sem subproject correspondente A campanha, NAO cai no primeiro da lista --
    # isso puxava o pipeline/metadados do SUB01 para dentro de outra campanha
    # (ex.: pipeline_SUB01.svg aparecendo no pacote do SUB04).
    subproject = next((s for s in subs if s.get("id") == campaign_id), {})
    active_sid = subproject.get("id")
    if include_pipeline and active_sid:
        r = client.get(f"/api/subproject/{active_sid}/pipeline.svg")
        if r.status_code == 200:
            (out / "assets").mkdir(parents=True, exist_ok=True)
            (out / "assets" / f"pipeline_{active_sid}.svg").write_bytes(r.get_data())

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
    config = {
        "title": (paper.get("config_title") if paper
                  else f"Docking Suite IFRJ · {campaign_id}"),
        "version": version,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "source": (paper.get("config_source") if paper
                   else "Resultados in silico (dados de trabalho)."),
        "compound": (paper.get("compound") if paper else None),
        "paper_title": subproject.get("titulo_paper"),
        "periodico_alvo": subproject.get("periodico_alvo"),
        "selectivity_cellular_max": (paper.get("selectivity_cellular_max") if paper else None),
        "ic50_note": (paper.get("ic50_note") if paper else None),
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


# Backward-compat: nome antigo. Prefira build_viewer(campaign_id).
def export_bundle(out_dir, compound="2-SCN-NQ", campaign_id="SUB01"):
    return build_viewer(campaign_id=campaign_id, out_dir=out_dir)


def main():
    ap = argparse.ArgumentParser(description="Empacota o NQ-Docking Viewer estatico (fiel ao paper congelado).")
    ap.add_argument("--out", default="dist_web")
    ap.add_argument("--campaign", default="SUB01", help="id da campanha (ex.: SUB01, SUB04)")
    ap.add_argument("--no-pipeline", action="store_true", help="nao gerar o SVG de pipeline")
    args = ap.parse_args()
    s = build_viewer(campaign_id=args.campaign, out_dir=args.out,
                     include_pipeline=not args.no_pipeline)
    tag = "fiel ao manuscrito submetido" if s["has_paper"] else "dados de trabalho (sem paper congelado)"
    print(f"[OK] bundle gerado ({tag}):")
    for k, v in s.items():
        print(f"     {k}: {v}")
    print(f"\nEnvie o conteudo de {s['out_dir']}/ para o Hostinger (ver README_DEPLOY.md).")


if __name__ == "__main__":
    main()
