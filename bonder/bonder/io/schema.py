#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""schema.py - contrato de dados unificado do bonder: o MSDataset.

Todo parser de motor de busca (MaxQuant/DIA-NN/FragPipe/generico) devolve UM
MSDataset com tres partes alinhadas pela mesma chave 'protein_id':

  intensity : DataFrame  (index = protein_id ; colunas = amostras ; valores float)
              AUSENTE = NaN, NUNCA 0. Chave = ACESSO UniProt lider (nunca gene).
  features  : DataFrame  (index = protein_id) com metadados por proteina:
              leading_accession, majority_protein_ids, gene_names, n_peptides,
              n_unique_peptides, q_value, is_contaminant, is_decoy, is_bait,
              species, shared_with_bait
  samples   : DataFrame  (index = sample) com o desenho experimental:
              condition, bait_id, role, replicate, raw_file

  provenance: dict com engine, engine_version, acquisition (DDA/DIA),
              quant_type (lfq/maxlfq/intensity), fdr_level, source_files, notes.

Serializa para Parquet (pandas + pyarrow) de forma reproduzivel.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas", "pyarrow"])
    import pandas as pd


# Colunas canonicas de features (ordem estavel).
FEATURE_COLUMNS = [
    "leading_accession",
    "majority_protein_ids",
    "gene_names",
    "n_peptides",
    "n_unique_peptides",
    "q_value",
    "is_contaminant",
    "is_decoy",
    "is_bait",
    "species",
    "shared_with_bait",
]

SAMPLE_COLUMNS = ["condition", "bait_id", "role", "replicate", "raw_file"]


def _parquet_available() -> bool:
    """True se houver engine de Parquet (pyarrow/fastparquet) REALMENTE importavel.

    Em Windows/conda o pyarrow as vezes instala mas falha no import (build/DLL);
    por isso testamos o import de verdade e caimos para pickle quando indisponivel,
    para a serializacao intermediaria nunca quebrar a analise.
    """
    for eng in ("pyarrow", "fastparquet"):
        try:
            __import__(eng)
            return True
        except Exception:
            continue
    return False


@dataclass
class MSDataset:
    """Matriz proteina x amostra unificada + metadados + proveniencia."""

    intensity: "pd.DataFrame"
    features: "pd.DataFrame"
    samples: "pd.DataFrame"
    provenance: Dict[str, Any] = field(default_factory=dict)

    # ---- validacao ----------------------------------------------------------
    def validate(self) -> "MSDataset":
        """Checa consistencia estrutural minima. Levanta ValueError se invalido."""
        if self.intensity.index.has_duplicates:
            dups = self.intensity.index[self.intensity.index.duplicated()].tolist()[:5]
            raise ValueError("protein_id duplicado na matriz de intensidade: %s ..." % dups)
        missing_feat = [c for c in FEATURE_COLUMNS if c not in self.features.columns]
        if missing_feat:
            raise ValueError("features sem colunas canonicas: %s" % missing_feat)
        # index de features deve cobrir a matriz
        not_covered = set(self.intensity.index) - set(self.features.index)
        if not_covered:
            raise ValueError("proteinas sem features: %s ..." % list(not_covered)[:5])
        # amostras devem bater com colunas da matriz
        cols = set(self.intensity.columns)
        smp = set(self.samples.index)
        if cols != smp:
            raise ValueError(
                "colunas da matriz e samples divergem. "
                "so na matriz=%s ; so em samples=%s"
                % (sorted(cols - smp)[:5], sorted(smp - cols)[:5])
            )
        return self

    # ---- consultas uteis para o scoring ------------------------------------
    def samples_in(self, condition: str) -> List[str]:
        s = self.samples
        return s.index[s["condition"] == condition].tolist()

    def observed_counts(self, condition: str) -> "pd.Series":
        """Numero de replicas com valor observado (nao-NaN) por proteina no grupo."""
        cols = self.samples_in(condition)
        if not cols:
            return pd.Series(0, index=self.intensity.index, dtype=int)
        return self.intensity[cols].notna().sum(axis=1).astype(int)

    def baits(self) -> List[str]:
        f = self.features
        return f.index[f["is_bait"].fillna(False).astype(bool)].tolist()

    def drop_contaminants_decoys(self) -> "MSDataset":
        """Remove decoys, contaminantes e 'only identified by site', preservando iscas."""
        f = self.features
        keep = ~(
            f["is_decoy"].fillna(False).astype(bool)
            | (f["is_contaminant"].fillna(False).astype(bool) & ~f["is_bait"].fillna(False).astype(bool))
        )
        idx = f.index[keep]
        return MSDataset(
            intensity=self.intensity.loc[idx].copy(),
            features=self.features.loc[idx].copy(),
            samples=self.samples.copy(),
            provenance={**self.provenance, "contaminants_removed": True},
        )

    # ---- (de)serializacao ---------------------------------------------------
    def to_parquet(self, run_dir: str | Path) -> Path:
        """Grava intensity/features/samples em Parquet + proveniencia em JSON."""
        d = Path(run_dir)
        d.mkdir(parents=True, exist_ok=True)
        fmt = "parquet" if _parquet_available() else "pickle"
        if fmt == "parquet":
            self.intensity.to_parquet(d / "intensity.parquet")
            self.features.to_parquet(d / "features.parquet")
            self.samples.to_parquet(d / "samples.parquet")
        else:
            # fallback resiliente quando pyarrow/fastparquet nao importam
            self.intensity.to_pickle(d / "intensity.pkl")
            self.features.to_pickle(d / "features.pkl")
            self.samples.to_pickle(d / "samples.pkl")
        prov = dict(self.provenance)
        prov["_storage_format"] = fmt
        (d / "provenance.json").write_text(
            json.dumps(prov, indent=2, ensure_ascii=True), encoding="utf-8"
        )
        return d

    @classmethod
    def from_parquet(cls, run_dir: str | Path) -> "MSDataset":
        d = Path(run_dir)
        prov = {}
        pj = d / "provenance.json"
        if pj.exists():
            prov = json.loads(pj.read_text(encoding="utf-8"))
        # detecta o formato pelo arquivo presente (parquet ou pickle)
        if (d / "intensity.parquet").exists():
            rd, ext = pd.read_parquet, "parquet"
        else:
            rd, ext = pd.read_pickle, "pkl"
        return cls(
            intensity=rd(d / ("intensity.%s" % ext)),
            features=rd(d / ("features.%s" % ext)),
            samples=rd(d / ("samples.%s" % ext)),
            provenance=prov,
        )

    def summary(self) -> Dict[str, Any]:
        f = self.features
        return {
            "n_proteins": int(self.intensity.shape[0]),
            "n_samples": int(self.intensity.shape[1]),
            "n_contaminants": int(f["is_contaminant"].fillna(False).astype(bool).sum()),
            "n_decoys": int(f["is_decoy"].fillna(False).astype(bool).sum()),
            "n_baits": int(f["is_bait"].fillna(False).astype(bool).sum()),
            "acquisition": self.provenance.get("acquisition"),
            "engine": self.provenance.get("engine"),
            "quant_type": self.provenance.get("quant_type"),
        }


def empty_features(index) -> "pd.DataFrame":
    """DataFrame de features vazio com as colunas canonicas e defaults sensatos."""
    df = pd.DataFrame(index=pd.Index(index, name="protein_id"), columns=FEATURE_COLUMNS)
    df["n_peptides"] = 0
    df["n_unique_peptides"] = 0
    df["q_value"] = pd.NA
    for c in ("is_contaminant", "is_decoy", "is_bait", "shared_with_bait"):
        df[c] = False
    df["species"] = "mouse"
    return df
