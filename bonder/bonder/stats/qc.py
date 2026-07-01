#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""qc.py - normalizacao e metricas de QC APENAS para visualizacao/dashboard.

IMPORTANTE: a inferencia estatistica (log2FC, p-valor, FDR) que vai para o paper
sai do DEP/limma/MSstats (R) e do SAINTexpress. As funcoes aqui servem para QC
exploratorio e graficos (log2, normalizacao por mediana, correlacao, contagem de
valores validos), NAO substituem a estatistica gold-standard.

Regra critica: NUNCA re-normalizar matrizes MaxLFQ do DIA-NN (ja normalizadas).
"""
from __future__ import annotations

import numpy as np

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["numpy", "pandas"])
    import pandas as pd


def log2_transform(intensity: "pd.DataFrame") -> "pd.DataFrame":
    return np.log2(intensity.where(intensity > 0))


def median_normalize(log_intensity: "pd.DataFrame", quant_type: str = "lfq") -> "pd.DataFrame":
    """Centraliza cada amostra na mediana global (so p/ QC). No-op para MaxLFQ."""
    if str(quant_type).lower() == "maxlfq":
        return log_intensity  # ja normalizado pelo DIA-NN
    col_med = log_intensity.median(axis=0, skipna=True)
    grand = col_med.median()
    return log_intensity.sub(col_med, axis=1).add(grand)


def valid_value_counts(ds) -> "pd.DataFrame":
    """Contagem de valores validos por amostra e por condicao (QC de completude)."""
    per_sample = ds.intensity.notna().sum(axis=0).rename("valid_values")
    df = per_sample.to_frame()
    df["condition"] = ds.samples["condition"].reindex(df.index).values
    return df


def sample_correlation(log_intensity: "pd.DataFrame") -> "pd.DataFrame":
    """Correlacao de Pearson entre amostras (QC de reprodutibilidade)."""
    return log_intensity.corr(method="pearson")


def qc_summary(ds) -> dict:
    log = log2_transform(ds.intensity)
    vv = valid_value_counts(ds)
    return {
        "n_proteins": int(ds.intensity.shape[0]),
        "n_samples": int(ds.intensity.shape[1]),
        "valid_values_per_sample": vv["valid_values"].to_dict(),
        "median_valid_values": float(vv["valid_values"].median()),
        "acquisition": ds.provenance.get("acquisition"),
        "quant_type": ds.provenance.get("quant_type"),
    }
