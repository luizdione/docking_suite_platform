#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""crapome.py - frequencia de background estilo CRAPome para AP-MS.

O CRAPome (Mellacheruvu et al., Nat Methods 2013) cataloga proteinas que aparecem
frequentemente em controles de AP-MS (background). Como a cobertura murina do
CRAPome publico e limitada, aqui priorizamos DUAS fontes:

  1. INTERNA: frequencia observada nos proprios controles do usuario
     (proteina vista em N de M corridas-controle) -> mais fiel ao seu sistema.
  2. EXTERNA: tabela CRAPome-style fornecida pelo usuario (opcional).

Preys com alta frequencia de background sao sinalizados (nao removidos): a decisao
final combina esse sinal com o escore SAINT/MSstats.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd


def internal_background(ds, control_conditions: Iterable[str]) -> "pd.DataFrame":
    """Frequencia de deteccao por proteina nos controles (observado em N de M)."""
    ctrl_cols = []
    for cond in control_conditions:
        ctrl_cols += ds.samples_in(cond)
    ctrl_cols = list(dict.fromkeys(ctrl_cols))  # unicos, ordem estavel
    m = len(ctrl_cols)
    if m == 0:
        return pd.DataFrame(
            {"n_observed": 0, "n_controls": 0, "background_freq": 0.0},
            index=ds.intensity.index,
        )
    observed = ds.intensity[ctrl_cols].notna().sum(axis=1).astype(int)
    df = pd.DataFrame({"n_observed": observed, "n_controls": m})
    df["background_freq"] = df["n_observed"] / m
    return df


def load_external_crapome(path, accession_col: str = "accession", freq_col: str = "background_freq") -> "pd.DataFrame":
    """Le uma tabela CRAPome-style do usuario (TSV) indexada por acesso."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if accession_col not in df.columns:
        raise ValueError("Tabela CRAPome sem coluna %r" % accession_col)
    df = df.set_index(accession_col)
    if freq_col in df.columns:
        df[freq_col] = pd.to_numeric(df[freq_col], errors="coerce")
    return df


def flag_background(scored: "pd.DataFrame", background: "pd.DataFrame", threshold: float = 0.5,
                    key: str = "prey") -> "pd.DataFrame":
    """Anexa background_freq e is_background aos resultados escoreados."""
    out = scored.copy()
    freq = background["background_freq"] if "background_freq" in background.columns else background.iloc[:, -1]
    out["background_freq"] = out[key].map(freq).fillna(0.0)
    out["is_background"] = out["background_freq"] >= threshold
    return out
