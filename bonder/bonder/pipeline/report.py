#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""report.py - montagem da tabela de interatores, figuras e saidas para a web.

A tabela final combina:
  - descritivos (media de intensidade e contagem de observados por grupo) - QC;
  - escore SAINT (saint_score, bfdr) - especificidade, quando disponivel;
  - DEP/limma (log2fc, p, FDR) - quantitativo, quando disponivel;
  - frequencia de background (controles) - filtro de contaminante.

Descritivos NAO substituem a estatistica: o ranqueamento prioriza BFDR (SAINT) e
FDR (DEP); os descritivos so entram como desempate/QC.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas", "matplotlib"])
    import pandas as pd


def assemble_interactors(
    ds,
    contrast: Dict[str, Any],
    saint_df: Optional["pd.DataFrame"] = None,
    dep_df: Optional["pd.DataFrame"] = None,
    background: Optional["pd.DataFrame"] = None,
    exclude_baits: bool = True,
) -> "pd.DataFrame":
    feats = ds.features
    idx = list(ds.intensity.index)
    if exclude_baits:
        bait_rows = set(feats.index[feats["is_bait"].fillna(False).astype(bool)])
        idx = [p for p in idx if p not in bait_rows]

    test_cols = ds.samples_in(contrast["test"])
    ctrl_cols = ds.samples_in(contrast["control"])
    mat = ds.intensity.loc[idx]
    mean_ip = mat[test_cols].mean(axis=1, skipna=True)
    mean_ctrl = mat[ctrl_cols].mean(axis=1, skipna=True)
    obs_ip = mat[test_cols].notna().sum(axis=1)
    obs_ctrl = mat[ctrl_cols].notna().sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2_fc_desc = np.log2(mean_ip / mean_ctrl)
    log2_fc_desc = log2_fc_desc.replace([np.inf, -np.inf], np.nan)

    df = pd.DataFrame({
        "prey": idx,
        "gene": feats.loc[idx, "gene_names"].astype(str).str.split(";").str[0].values,
        "mean_intensity_ip": mean_ip.values,
        "mean_intensity_ctrl": mean_ctrl.values,
        "observed_ip": obs_ip.values,
        "observed_ctrl": obs_ctrl.values,
        "log2fc_descriptive": log2_fc_desc.values,
    })

    if saint_df is not None and not saint_df.empty:
        s = saint_df[[c for c in ("prey", "saint_score", "bfdr", "fold_change") if c in saint_df.columns]]
        df = df.merge(s, on="prey", how="left")
    if dep_df is not None and not dep_df.empty:
        d = dep_df.rename(columns={"accession": "prey"})
        keep = [c for c in ("prey", "log2fc", "p_value", "adj_p_value", "significant") if c in d.columns]
        df = df.merge(d[keep], on="prey", how="left")
    if background is not None and not background.empty and "background_freq" in background.columns:
        df["background_freq"] = df["prey"].map(background["background_freq"]).fillna(0.0)

    # ranqueamento: BFDR (SAINT) -> FDR (DEP) -> fold change descritivo
    sort_cols, ascending = [], []
    if "bfdr" in df.columns:
        sort_cols.append("bfdr"); ascending.append(True)
    if "adj_p_value" in df.columns:
        sort_cols.append("adj_p_value"); ascending.append(True)
    sort_cols.append("log2fc_descriptive"); ascending.append(False)
    df = df.sort_values(sort_cols, ascending=ascending, na_position="last").reset_index(drop=True)

    df["confidence"] = _confidence_label(df)
    return df


def _confidence_label(df: "pd.DataFrame") -> "pd.Series":
    label = pd.Series("baixa", index=df.index)
    if "bfdr" in df.columns:
        label = label.mask(df["bfdr"] <= 0.05, "alta")
        label = label.mask((df["bfdr"] > 0.05) & (df["bfdr"] <= 0.10), "media")
    elif "adj_p_value" in df.columns:
        lfc = df["log2fc"].abs() if "log2fc" in df.columns else pd.Series(0.0, index=df.index)
        label = label.mask((df["adj_p_value"] <= 0.05) & (lfc >= 1), "alta")
    return label


def volcano_plot(dep_df: "pd.DataFrame", out_path, alpha: float = 0.05, lfc: float = 1.0) -> Optional[Path]:
    """Volcano log2FC x -log10(FDR) a partir da saida DEP/limma."""
    if dep_df is None or dep_df.empty or "log2fc" not in dep_df.columns:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = pd.to_numeric(dep_df["log2fc"], errors="coerce")
    p = pd.to_numeric(dep_df.get("adj_p_value"), errors="coerce")
    y = -np.log10(p.replace(0, np.nan))
    sig = (p <= alpha) & (x.abs() >= lfc)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x[~sig], y[~sig], s=8, c="#9aa0a6", alpha=0.6, label="ns")
    ax.scatter(x[sig], y[sig], s=12, c="#d93025", alpha=0.8, label="enriquecido")
    ax.axhline(-np.log10(alpha), ls="--", c="#5f6368", lw=0.8)
    ax.axvline(lfc, ls="--", c="#5f6368", lw=0.8)
    ax.axvline(-lfc, ls="--", c="#5f6368", lw=0.8)
    ax.set_xlabel("log2 fold-change (IP / controle)")
    ax.set_ylabel("-log10 FDR")
    ax.set_title("Enriquecimento AP-MS (DEP/limma)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def _sanitize(obj):
    """Converte floats nao-finitos (NaN/inf) em None para JSON valido no navegador."""
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else None
    # numpy scalars
    if isinstance(obj, np.floating):
        f = float(obj)
        return f if np.isfinite(f) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def write_results_json(run_dir, payload: Dict[str, Any]) -> Path:
    import json
    p = Path(run_dir) / "results.json"
    clean = _sanitize(payload)
    # allow_nan=False garante JSON estrito (JSON.parse do navegador rejeita NaN)
    p.write_text(json.dumps(clean, indent=2, ensure_ascii=True, default=str, allow_nan=False),
                 encoding="utf-8")
    return p
