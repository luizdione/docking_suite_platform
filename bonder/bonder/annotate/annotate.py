#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""annotate.py - anotacao dos interatores e enriquecimento funcional.

Anota o resultado com gene/nome (UniProt murino), marca parceiros ja conhecidos
(STRING/BioGRID) e, opcionalmente, roda enriquecimento GO/via com g:Profiler.
O enriquecimento e opcional e faz chamada de rede; sem a lib ou offline, retorna
None com aviso (nao quebra o pipeline).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd

from bonder.db import string_biogrid


def annotate_results(
    results: "pd.DataFrame",
    mouse_ref=None,
    edges: Optional["pd.DataFrame"] = None,
    bait_accession: Optional[str] = None,
    accession_col: str = "prey",
) -> "pd.DataFrame":
    """Adiciona gene/known_interactor aos resultados escoreados."""
    out = results.copy()
    accs = out[accession_col].astype(str)
    if mouse_ref is not None:
        gene_map = mouse_ref.gene_map()
        if gene_map:
            out["gene"] = accs.map(lambda a: gene_map.get(a)).fillna(out.get("prey_gene"))
    if edges is not None and bait_accession is not None:
        known = string_biogrid.known_interactors(edges, bait_accession)
        out["known_interactor"] = accs.isin(known)
    return out


def go_enrichment(genes: List[str], organism: str = "mmusculus") -> Optional["pd.DataFrame"]:
    """Enriquecimento GO/via via g:Profiler (opcional). None se indisponivel."""
    try:
        from gprofiler import GProfiler
    except ImportError:
        return None
    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.profile(organism=organism, query=[g for g in genes if g])
        return res
    except Exception:
        return None
