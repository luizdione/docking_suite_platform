#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Integracao ponta-a-ponta do pipeline em CPU (SAINT/R ausentes -> pulados)."""
import json
from pathlib import Path

import pandas as pd

from bonder.config import load_config, PROJECT_ROOT


def _write_maxquant_full(path):
    # colunas casam com samples.yaml (NSs_IP_1..3, NSs_CTRL_1..3)
    cols = {
        "Majority protein IDs": ["P16110", "NSs_ICOARACI", "P12345", "Q8R1", "CON__P00761"],
        "Gene names": ["Lgals3", "NSs", "Prey1", "Prey2", "Trypsin"],
        "Peptides": ["8", "12", "5", "6", "3"],
        "Razor + unique peptides": ["6", "10", "4", "5", "2"],
        "Reverse": ["", "", "", "", ""],
        "Potential contaminant": ["", "", "", "", "+"],
    }
    for s, base in [("NSs_IP_1", 5000), ("NSs_IP_2", 5200), ("NSs_IP_3", 4800)]:
        cols["LFQ intensity %s" % s] = [str(base + 100), str(base * 10), str(base), str(base // 2), "300"]
    for s in ["NSs_CTRL_1", "NSs_CTRL_2", "NSs_CTRL_3"]:
        cols["LFQ intensity %s" % s] = ["0", "50", "700", "0", "320"]
    pd.DataFrame(cols).to_csv(path, sep="\t", index=False)


def test_pipeline_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setenv("BONDER_CONFIG_DIR", str(PROJECT_ROOT / "config"))
    monkeypatch.setenv("BONDER_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("BONDER_DATA_ROOT", str(tmp_path / "data"))
    cfg = load_config()

    mq = tmp_path / "proteinGroups.txt"
    _write_maxquant_full(mq)

    from bonder.pipeline import run_pipeline
    contrast = next(c for c in cfg.contrasts() if c["name"] == "NSs_vs_CTRL")
    run_dir = run_pipeline(mq, contrast, "test_run", cfg=cfg)

    run_dir = Path(run_dir)
    for artifact in ("interactors.tsv", "manifest.json", "qc.json", "network.json", "results.json"):
        assert (run_dir / artifact).exists(), "faltou %s" % artifact
    # SAINT/DEP inputs sempre gerados mesmo sem os binarios
    assert (run_dir / "saint" / "interaction.txt").exists()
    assert (run_dir / "dep" / "data.tsv").exists()

    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["provenance"]["engine"] == "maxquant"
    assert manifest["device"]["device"] == "cpu"  # sem GPU no CI
    steps = {s["step"]: s["status"] for s in manifest["steps"]}
    assert steps["saint"] in ("pulado", "erro")  # binario ausente
    assert steps["dep"] in ("pulado", "erro")     # Rscript ausente
    assert steps["ingest"] == "ok"

    table = pd.read_csv(run_dir / "interactors.tsv", sep="\t")
    # isca exogena e contaminante nao entram como candidatos
    assert "NSs_ICOARACI" not in set(table["prey"])
    assert "P00761" not in set(table["prey"])
    # prey real presente
    assert "P12345" in set(table["prey"])
