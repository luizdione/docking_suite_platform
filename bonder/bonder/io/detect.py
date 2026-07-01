#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""detect.py - deteccao automatica de motor de busca e modo de aquisicao.

O usuario NAO precisa dizer se a Co-IP foi DDA ou DIA: inferimos do proprio
arquivo (nome + cabecalho de colunas) e roteamos para o parser correto.

  MaxQuant  proteinGroups.txt      -> DDA
  DIA-NN    report.pg_matrix.tsv   -> DIA
  FragPipe  combined_protein.tsv   -> DDA (ou DIA se detectado)
  generico  requer sidecar YAML explicito

Retorna (engine, acquisition). load_search_output() ja devolve o MSDataset.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd

from bonder.io import parsers


def _header(path: Path) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.readline()


def detect_engine(path) -> Tuple[str, str]:
    """Infere (engine, acquisition) a partir do nome do arquivo e do cabecalho."""
    path = Path(path)
    name = path.name.lower()

    # 1) por nome de arquivo (sinal forte)
    if "proteingroups" in name:
        return "maxquant", "DDA"
    if "pg_matrix" in name or ("report" in name and "diann" in name):
        return "diann", "DIA"
    if "combined_protein" in name:
        return "fragpipe", "DDA"

    # 2) por cabecalho de colunas
    try:
        header = _header(path)
    except OSError:
        header = ""
    cols = {c.strip() for c in header.rstrip("\n").split("\t")}

    if {"Majority protein IDs", "Reverse"} & cols or "Potential contaminant" in cols:
        return "maxquant", "DDA"
    if "Protein.Group" in cols or "Protein.Ids" in cols:
        return "diann", "DIA"
    if "Protein ID" in cols and any(c.endswith("Intensity") for c in cols):
        return "fragpipe", "DDA"

    raise ValueError(
        "Motor de busca nao reconhecido para %s. Use um sidecar generico "
        "(parse_generic) descrevendo id_col e as colunas de intensidade." % path
    )


def load_search_output(
    path,
    samples_cfg=None,
    baits_cfg=None,
    sidecar: Optional[Dict[str, Any]] = None,
    report_tsv=None,
) -> "parsers.MSDataset":
    """Ponto de entrada unico: detecta o motor e devolve o MSDataset unificado."""
    if sidecar is not None:
        return parsers.parse_generic(path, sidecar, samples_cfg=samples_cfg, baits_cfg=baits_cfg)

    engine, acquisition = detect_engine(path)
    if engine == "maxquant":
        return parsers.parse_maxquant(path, samples_cfg=samples_cfg, baits_cfg=baits_cfg)
    if engine == "diann":
        return parsers.parse_diann(
            path, report_tsv=report_tsv, samples_cfg=samples_cfg, baits_cfg=baits_cfg
        )
    if engine == "fragpipe":
        return parsers.parse_fragpipe(
            path, samples_cfg=samples_cfg, baits_cfg=baits_cfg, acquisition=acquisition
        )
    raise ValueError("engine detectado sem parser: %s" % engine)
