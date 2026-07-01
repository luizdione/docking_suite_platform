#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""uniprot.py - proteoma murino UniProt + crosswalks acesso<->gene/STRING/BioGRID.

Le o FASTA (sequencias, comprimentos) e a tabela TSV de id-mapping do UniProtKB
(taxid 10090). Fornece um objeto MouseReference com consultas usadas pelo scoring
(comprimento de prey para o SAINT), anotacao (gene/nome) e rede (STRING/BioGRID).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd


def iter_fasta(path) -> Iterator[Tuple[str, str, str]]:
    """Itera (accession, header, sequence) de um FASTA UniProt (sp|ACC|NAME ...)."""
    acc = header = None
    seq: list = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    yield acc, header, "".join(seq)
                header = line[1:]
                parts = header.split("|")
                acc = parts[1] if len(parts) >= 2 else header.split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if acc is not None:
            yield acc, header, "".join(seq)


def read_fasta(path) -> Dict[str, str]:
    """accession -> sequencia."""
    return {acc: seq for acc, _h, seq in iter_fasta(path)}


def sequence_lengths(path) -> Dict[str, int]:
    """accession -> comprimento da sequencia (residuos)."""
    return {acc: len(seq) for acc, _h, seq in iter_fasta(path)}


class MouseReference:
    """Consultas sobre o proteoma murino de referencia."""

    def __init__(
        self,
        lengths: Optional[Dict[str, int]] = None,
        idmap: Optional["pd.DataFrame"] = None,
    ):
        self.lengths = lengths or {}
        # idmap indexado por accession
        self.idmap = idmap if idmap is not None else pd.DataFrame()

    @classmethod
    def from_files(cls, fasta_path=None, tsv_path=None) -> "MouseReference":
        lengths = sequence_lengths(fasta_path) if fasta_path and Path(fasta_path).exists() else {}
        idmap = pd.DataFrame()
        if tsv_path and Path(tsv_path).exists():
            idmap = pd.read_csv(tsv_path, sep="\t", dtype=str)
            # UniProt REST usa 'Entry' como acesso quando fields=accession -> normaliza
            acc_col = "Entry" if "Entry" in idmap.columns else idmap.columns[0]
            idmap = idmap.rename(columns={acc_col: "accession"}).set_index("accession")
        return cls(lengths=lengths, idmap=idmap)

    def length_for(self, accession: str, default: int = 500) -> int:
        return int(self.lengths.get(accession, default))

    def gene_for(self, accession: str) -> Optional[str]:
        if accession in self.idmap.index:
            for col in ("Gene Names (primary)", "gene_primary", "Gene names (primary)"):
                if col in self.idmap.columns:
                    val = self.idmap.at[accession, col]
                    if isinstance(val, str) and val.strip():
                        return val.split(";")[0].strip()
        return None

    def gene_map(self) -> Dict[str, str]:
        """Mapa accession -> gene primario, para uso em massa (SAINT/anotacao)."""
        out: Dict[str, str] = {}
        if self.idmap.empty:
            return out
        for col in ("Gene Names (primary)", "gene_primary", "Gene names (primary)"):
            if col in self.idmap.columns:
                for acc, val in self.idmap[col].items():
                    if isinstance(val, str) and val.strip():
                        out[acc] = val.split(";")[0].strip()
                break
        return out
