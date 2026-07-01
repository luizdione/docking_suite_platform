#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""orchestrator.py - unico coordenador cross-layer do bonder.

Fluxo por contraste (isca vs controle):
  1. ingesta        io.load_search_output -> MSDataset (detecta motor/aquisicao)
  2. QC             stats.qc.qc_summary
  3. limpeza        remove decoys/contaminantes (preserva iscas)
  4. especificidade scoring.SAINT (e MiST p/ isca viral) - se binario disponivel
  5. quantitativo   stats DEP/limma (R) - se Rscript disponivel
  6. background     CRAPome interno (controles do usuario)
  7. montagem       report.assemble_interactors (ranqueia por BFDR/FDR)
  8. anotacao+rede  annotate + network
  9. saidas         TSV, volcano, graphml/cytoscape json, results.json
 10. manifesto      manifest.json (versoes, ferramentas, device, checksums)

Etapas cujas ferramentas externas faltam (SAINT/R/GPU) sao PULADAS com registro
explicito no manifesto - o pipeline nunca finge cobertura que nao teve.
"""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from bonder import __version__
from bonder.config import BonderConfig, load_config
from bonder.io import load_search_output
from bonder.stats import qc
from bonder.db import crapome
from bonder.db.uniprot import MouseReference
from bonder.pipeline import report


def _tool_available(name: str) -> bool:
    import shutil
    return shutil.which(name) is not None


def run_pipeline(
    input_path,
    contrast: Dict[str, Any],
    run_id: str,
    cfg: Optional[BonderConfig] = None,
    report_tsv=None,
    enable_gpu_esm2: bool = False,
) -> Path:
    """Executa o pipeline para um contraste. Retorna o diretorio da corrida."""
    cfg = cfg or load_config()
    run_dir = cfg.run_path(run_id)
    steps: List[Dict[str, Any]] = []

    def log_step(name, status, **extra):
        steps.append({"step": name, "status": status, **extra})

    # 1) ingesta ------------------------------------------------------------
    baits_cfg = cfg.baits
    samples_cfg = cfg.samples.get("samples", [])
    ds = load_search_output(input_path, samples_cfg=samples_cfg, baits_cfg=baits_cfg,
                            report_tsv=report_tsv)
    ds.to_parquet(run_dir / "dataset")
    log_step("ingest", "ok", **ds.summary())

    # 2) QC -----------------------------------------------------------------
    qc_sum = qc.qc_summary(ds)
    (run_dir / "qc.json").write_text(json.dumps(qc_sum, indent=2, ensure_ascii=True), encoding="utf-8")
    log_step("qc", "ok", median_valid_values=qc_sum["median_valid_values"])

    # 3) limpeza ------------------------------------------------------------
    clean = ds.drop_contaminants_decoys()
    log_step("clean", "ok", n_proteins=int(clean.intensity.shape[0]))

    # referencia murina (opcional; usada p/ comprimento de prey e genes)
    mouse_ref = None
    fasta = cfg.cache_path("uniprot", "mouse.fasta")
    tsv = cfg.cache_path("uniprot", "mouse_idmap.tsv")
    if Path(fasta).exists() or Path(tsv).exists():
        mouse_ref = MouseReference.from_files(
            fasta if Path(fasta).exists() else None,
            tsv if Path(tsv).exists() else None,
        )

    # 4) especificidade: SAINT (+ MiST p/ isca viral) -----------------------
    saint_df = None
    saint_dir = run_dir / "saint"
    from bonder.scoring import saint as saint_mod, mist as mist_mod
    prey_lengths = mouse_ref.lengths if mouse_ref else None
    gene_map = mouse_ref.gene_map() if mouse_ref else None
    inputs = saint_mod.write_saint_inputs(clean, contrast, saint_dir,
                                          prey_lengths=prey_lengths, gene_map=gene_map)
    if _tool_available("SAINTexpress-int"):
        try:
            lst = saint_mod.run_saintexpress(inputs, saint_dir)
            saint_df = saint_mod.parse_saint_output(lst)
            log_step("saint", "ok", n_prey=int(saint_df.shape[0]))
        except Exception as e:
            log_step("saint", "erro", detail=str(e))
    else:
        log_step("saint", "pulado", detail="SAINTexpress-int ausente; inputs gerados em saint/")

    bait = cfg.bait(contrast["bait_id"]) if contrast.get("bait_id") in cfg.bait_ids() else {}
    if (bait.get("scoring") or {}).get("mist"):
        mist_inputs = mist_mod.write_mist_inputs(clean, contrast, run_dir / "mist")
        mist_out = None
        try:
            mist_out = mist_mod.run_mist(mist_inputs, run_dir / "mist")
        except Exception as e:
            log_step("mist", "erro", detail=str(e))
        log_step("mist", "ok" if mist_out else "pulado",
                 detail="MiST executado" if mist_out else "BONDER_MIST_CMD nao definido; inputs gerados")

    # 5) quantitativo: DEP/limma (R) ---------------------------------------
    dep_df = None
    from bonder.stats import prepare, runner
    dep_dir = run_dir / "dep"
    dep_inputs = prepare.write_dep_inputs(clean, contrast, dep_dir)
    if _tool_available("Rscript"):
        try:
            out_tsv = runner.run_dep(dep_inputs, dep_dir)
            dep_df = runner.parse_dep_results(out_tsv, contrast["test"], contrast["control"])
            log_step("dep", "ok", n_significant=int(dep_df.get("significant", []).sum()) if dep_df is not None else 0)
        except Exception as e:
            log_step("dep", "erro", detail=str(e))
    else:
        log_step("dep", "pulado", detail="Rscript ausente; inputs DEP gerados em dep/")

    # 6) background (CRAPome interno) --------------------------------------
    bg = crapome.internal_background(ds, [contrast["control"]])
    log_step("background", "ok")

    # 7) montagem ----------------------------------------------------------
    table = report.assemble_interactors(clean, contrast, saint_df=saint_df, dep_df=dep_df, background=bg)
    table_path = run_dir / "interactors.tsv"
    table.to_csv(table_path, sep="\t", index=False)
    log_step("assemble", "ok", n_candidates=int(table.shape[0]))

    # 8) anotacao + rede ---------------------------------------------------
    from bonder.annotate import build_interactome, to_cytoscape_json, export_graphml
    top = table.head(150)  # rede centrada nos melhores candidatos
    g = build_interactome(contrast["bait_id"], top)
    export_graphml(g, run_dir / "network.graphml")
    (run_dir / "network.json").write_text(
        json.dumps(to_cytoscape_json(g), indent=2, ensure_ascii=True), encoding="utf-8"
    )
    log_step("network", "ok", n_nodes=g.number_of_nodes(), n_edges=g.number_of_edges())

    # 9) figuras + results.json -------------------------------------------
    volcano = report.volcano_plot(dep_df, run_dir / "volcano.png")
    payload = {
        "run_id": run_id,
        "contrast": contrast,
        "summary": ds.summary(),
        "qc": qc_sum,
        "steps": steps,
        "top_interactors": table.head(50).to_dict(orient="records"),
        "artifacts": {
            "interactors_tsv": str(table_path.name),
            "network_graphml": "network.graphml",
            "network_json": "network.json",
            "volcano_png": volcano.name if volcano else None,
        },
    }
    report.write_results_json(run_dir, payload)

    # 10) manifesto de reprodutibilidade ----------------------------------
    from bonder.gpu.detect import probe
    manifest = {
        "bonder_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "input": str(input_path),
        "contrast": contrast,
        "tools": {
            "SAINTexpress-int": _tool_available("SAINTexpress-int"),
            "Rscript": _tool_available("Rscript"),
            "diann": _tool_available("diann"),
            "colabfold_batch": _tool_available("colabfold_batch"),
        },
        "device": probe().to_dict(),
        "provenance": ds.provenance,
        "steps": steps,
    }
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    return run_dir
