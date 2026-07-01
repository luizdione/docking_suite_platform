#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""runner.py - execucao dos scripts R de estatistica (DEP/limma, MSstats).

Envolve o Rscript. Os scripts R vivem em bonder/R/ e rodam no conteiner de
estatistica (Dockerfile.r). Em ambiente sem R, run_dep levanta erro claro; a
preparacao de entradas (prepare.py) e a leitura de saidas continuam testaveis.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd

R_DIR = Path(__file__).resolve().parents[2] / "R"


def _rscript() -> str:
    exe = shutil.which("Rscript")
    if exe is None:
        raise FileNotFoundError(
            "Rscript nao encontrado. A estatistica gold-standard (DEP/limma/MSstats) "
            "roda no conteiner de estatistica: use 'docker compose run stats' ou "
            "instale R + Bioconductor (DEP, limma, MSstats)."
        )
    return exe


def run_dep(inputs: Dict[str, Any], out_dir, alpha: float = 0.05, lfc: float = 1.0) -> Path:
    """Executa run_dep_limma.R. Retorna o caminho do results.tsv."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "dep_results.tsv"
    cmd = [
        _rscript(), str(R_DIR / "run_dep_limma.R"),
        str(inputs["data"]), str(inputs["design"]), str(out_tsv),
        str(inputs["test"]), str(inputs["control"]),
        str(alpha), str(lfc),
    ]
    subprocess.run(cmd, check=True)
    if not out_tsv.exists():
        raise RuntimeError("run_dep_limma.R terminou sem gerar %s" % out_tsv)
    return out_tsv


def parse_dep_results(out_tsv, test: str, control: str) -> "pd.DataFrame":
    """Padroniza a saida do DEP: accession, gene, log2fc, p_value, adj_p_value, significant."""
    df = pd.read_csv(out_tsv, sep="\t", dtype=str)
    contrast = "%s_vs_%s" % (test, control)
    ratio_c = "%s_ratio" % contrast
    p_c = "%s_p.val" % contrast
    padj_c = "%s_p.adj" % contrast
    sig_c = "%s_significant" % contrast

    def _col(name, alt=None):
        if name in df.columns:
            return df[name]
        if alt and alt in df.columns:
            return df[alt]
        return pd.Series([None] * len(df))

    out = pd.DataFrame({
        "accession": _col("ID", "name"),
        "gene": _col("name"),
        "log2fc": pd.to_numeric(_col(ratio_c), errors="coerce"),
        "p_value": pd.to_numeric(_col(p_c), errors="coerce"),
        "adj_p_value": pd.to_numeric(_col(padj_c), errors="coerce"),
    })
    if sig_c in df.columns:
        out["significant"] = df[sig_c].astype(str).str.upper().isin(["TRUE", "1", "YES"])
    else:
        out["significant"] = (out["adj_p_value"] < 0.05) & (out["log2fc"].abs() >= 1)
    return out
