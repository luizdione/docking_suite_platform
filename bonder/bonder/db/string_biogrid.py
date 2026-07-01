#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""string_biogrid.py - carregadores das redes PPI murinas STRING e BioGRID.

STRING v12.0 (10090): 'protein.links' + 'protein.aliases' para mapear os ids
STRING (10090.ENSMUSP...) para acessos UniProt. Arestas com combined_score abaixo
do limiar sao descartadas. Arestas cujo id nao mapear ficam auditaveis (retornadas
com uniprot None em vez de descartadas silenciosamente).

BioGRID: arquivo tab3 do organismo, usando as colunas 'SWISS-PROT Accessions'.

Saida padrao: DataFrame com colunas [a, b, score, source].
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def string_alias_map(aliases_path) -> Dict[str, str]:
    """id STRING -> acesso UniProt (via linhas de alias com source UniProt)."""
    amap: Dict[str, str] = {}
    with _open_text(aliases_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            string_id, alias, source = parts[0], parts[1], parts[2]
            if "UniProt_AC" in source or "Uniprot_AC" in source:
                amap.setdefault(string_id, alias)
    return amap


def load_string_edges(links_path, aliases_path, min_score: int = 400) -> "pd.DataFrame":
    """Arestas STRING mapeadas a UniProt, filtradas por combined_score >= min_score."""
    amap = string_alias_map(aliases_path)
    rows = []
    with _open_text(links_path) as fh:
        header = fh.readline()  # 'protein1 protein2 combined_score'
        for line in fh:
            p1, p2, score = line.split()
            score = int(score)
            if score < min_score:
                continue
            rows.append((amap.get(p1), amap.get(p2), score, "STRING"))
    return pd.DataFrame(rows, columns=["a", "b", "score", "source"])


def load_biogrid_edges(tab3_path) -> "pd.DataFrame":
    """Arestas BioGRID (Mus musculus) por acesso Swiss-Prot."""
    df = pd.read_csv(tab3_path, sep="\t", dtype=str, low_memory=False)
    a_col = "SWISS-PROT Accessions Interactor A"
    b_col = "SWISS-PROT Accessions Interactor B"
    if a_col not in df.columns or b_col not in df.columns:
        raise ValueError("BioGRID tab3 sem colunas SWISS-PROT esperadas: %s" % tab3_path)
    out = df[[a_col, b_col]].copy()
    out.columns = ["a", "b"]
    out["a"] = out["a"].str.split("|").str[0]
    out["b"] = out["b"].str.split("|").str[0]
    out = out[(out["a"] != "-") & (out["b"] != "-")].dropna()
    out["score"] = 1000  # BioGRID = evidencia experimental curada
    out["source"] = "BioGRID"
    return out.reset_index(drop=True)


def known_interactors(edges: "pd.DataFrame", accession: str) -> set:
    """Conjunto de parceiros conhecidos de `accession` na rede fornecida."""
    if edges.empty:
        return set()
    a = set(edges.loc[edges["a"] == accession, "b"].dropna())
    b = set(edges.loc[edges["b"] == accession, "a"].dropna())
    return a | b
