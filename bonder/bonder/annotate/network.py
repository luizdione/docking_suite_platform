#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""network.py - montagem da rede de interacao centrada na isca.

Combina os interatores detectados (SAINT/DEP) com as arestas de referencia
(STRING/BioGRID) para desenhar um interactoma centrado na isca. Exporta GraphML
(Cytoscape) e JSON (cytoscape.js/dashboard).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import networkx as nx
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["networkx", "pandas"])
    import networkx as nx
    import pandas as pd


def build_interactome(
    bait_id: str,
    interactors: "pd.DataFrame",
    edges: Optional["pd.DataFrame"] = None,
    accession_col: str = "prey",
    score_cols: Optional[List[str]] = None,
) -> "nx.Graph":
    """Grafo isca-presas + arestas presa-presa conhecidas (STRING/BioGRID).

    interactors: DataFrame com uma coluna de acesso (accession_col) e colunas de
    escore (saint_score, bfdr, log2fc, background_freq, prey_gene...).
    edges: DataFrame [a, b, score, source] (opcional).
    """
    score_cols = score_cols or ["saint_score", "bfdr", "log2fc", "background_freq", "prey_gene", "gene"]
    g = nx.Graph()
    g.add_node(bait_id, node_type="bait", is_bait=True, label=bait_id)

    preys = set()
    for _i, row in interactors.iterrows():
        acc = str(row[accession_col])
        preys.add(acc)
        attrs: Dict[str, Any] = {"node_type": "prey", "is_bait": False, "label": acc}
        for c in score_cols:
            if c in interactors.columns and pd.notna(row.get(c)):
                val = row[c]
                attrs[c] = float(val) if isinstance(val, (int, float)) else str(val)
        if attrs.get("gene") or attrs.get("prey_gene"):
            attrs["label"] = attrs.get("prey_gene") or attrs.get("gene")
        g.add_node(acc, **attrs)
        g.add_edge(bait_id, acc, edge_type="detected", source="bonder")

    # arestas presa-presa a partir da referencia (apenas entre preys detectados)
    if edges is not None and not edges.empty:
        for _i, e in edges.iterrows():
            a, b = e.get("a"), e.get("b")
            if a in preys and b in preys and a != b:
                g.add_edge(a, b, edge_type="reference",
                           source=str(e.get("source", "STRING")),
                           score=float(e.get("score", 0) or 0))
    return g


def to_cytoscape_json(g: "nx.Graph") -> Dict[str, Any]:
    """Serializa o grafo para o formato elements do cytoscape.js."""
    nodes = [{"data": {"id": str(n), **{k: v for k, v in d.items()}}} for n, d in g.nodes(data=True)]
    edges = []
    for u, v, d in g.edges(data=True):
        edges.append({"data": {"source": str(u), "target": str(v), **{k: val for k, val in d.items()}}})
    return {"elements": {"nodes": nodes, "edges": edges}}


def export_graphml(g: "nx.Graph", path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(g, path)
    return path
