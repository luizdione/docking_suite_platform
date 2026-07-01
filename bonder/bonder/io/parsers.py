#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""parsers.py - parsers por motor de busca -> MSDataset unificado.

Um parser por motor, com MAPA DE COLUNAS EXPLICITO (sem autodeteccao de coluna).
Regras invariantes (ver bonder/io/schema.py):
  - chave = ACESSO UniProt lider (primeiro token do grupo de proteinas)
  - ausente = NaN, nunca 0 (0 do MaxQuant/FragPipe -> NaN em DDA)
  - decoys/contaminantes flaggeados por coluna booleana E por prefixo de acesso
  - iscas EXOGENAS marcadas (is_bait) e mantidas distintas do ortologo murino

parse_maxquant / parse_diann / parse_fragpipe / parse_generic recebem:
  path         : arquivo de saida do motor
  samples_cfg  : lista de dicts de samples.yaml (para reconciliar colunas) ou None
  baits_cfg    : conteudo de baits.yaml (para marcar iscas) ou None
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["numpy", "pandas", "pyarrow"])
    import numpy as np
    import pandas as pd

from bonder.io import contaminants as cont
from bonder.io.schema import MSDataset, empty_features


# ---------------------------------------------------------------------------
# Helpers compartilhados
# ---------------------------------------------------------------------------
def _leading_accession(protein_group: Any) -> str:
    """Primeiro token de um grupo de proteinas (separado por ';'), sem prefixos."""
    if not isinstance(protein_group, str) or not protein_group.strip():
        return ""
    first = protein_group.split(";")[0].strip()
    # formato 'sp|P17931|LEG3_HUMAN' -> pega o acesso do meio quando presente
    if "|" in first:
        parts = first.split("|")
        if len(parts) >= 2 and parts[1]:
            first = parts[1]
    return cont.strip_prefixes(first)


def _to_num(series, zero_to_nan: bool) -> "np.ndarray":
    """Converte para float64 posicional (numpy). Retorna array para evitar
    realinhamento por index ao atribuir a colunas de DataFrame."""
    s = pd.to_numeric(pd.Series(series), errors="coerce")
    if zero_to_nan:
        s = s.mask(s == 0.0, other=np.nan)
    return s.astype("float64").to_numpy()


def _num_col(df: "pd.DataFrame", names, n: int, default: int = 0) -> "np.ndarray":
    """Primeira coluna numerica existente dentre `names`; senao vetor de default."""
    if isinstance(names, str):
        names = [names]
    for nm in names:
        if nm in df.columns:
            return pd.to_numeric(df[nm], errors="coerce").fillna(default).astype(int).values
    return np.full(n, default, dtype=int)


def _str_col(df: "pd.DataFrame", names, n: int, default: str = "") -> "np.ndarray":
    if isinstance(names, str):
        names = [names]
    for nm in names:
        if nm in df.columns:
            return df[nm].fillna(default).astype(str).values
    return np.full(n, default, dtype=object)


def _match_intensity_columns(columns: List[str], regexes: List[str]) -> Tuple[List[str], List[str], str]:
    """Casa colunas de intensidade pela primeira regex que der match.

    Retorna (colunas_originais, rotulos_run_capturados, regex_usada).
    """
    for rx in regexes:
        pat = re.compile(rx)
        matched = [(c, pat.match(c)) for c in columns]
        hits = [(c, m.group(1).strip()) for c, m in matched if m]
        if hits:
            orig = [c for c, _ in hits]
            runs = [r for _, r in hits]
            return orig, runs, rx
    return [], [], ""


def _reconcile_samples(
    run_labels: List[str], samples_cfg: Optional[List[Dict[str, Any]]]
) -> Tuple[Dict[str, str], "pd.DataFrame"]:
    """Mapeia rotulos de run detectados -> sample ids e monta o samples DataFrame.

    Estrategia de casamento (na ordem): sample == run ; stem(raw_file) == run ;
    raw_file == run ; run contido em raw_file/sample. Erra se algum run ficar sem
    mapeamento (evita atribuir replica ao grupo errado silenciosamente).
    """
    if not samples_cfg:
        # Sem metadados: usa o proprio rotulo como sample e condicao.
        rename = {r: r for r in run_labels}
        df = pd.DataFrame(
            {
                "condition": run_labels,
                "bait_id": [None] * len(run_labels),
                "role": [None] * len(run_labels),
                "replicate": list(range(1, len(run_labels) + 1)),
                "raw_file": run_labels,
            },
            index=pd.Index(run_labels, name="sample"),
        )
        return rename, df

    def _stem(x: Any) -> str:
        return Path(str(x)).stem if x is not None else ""

    rename: Dict[str, str] = {}
    used_samples: set = set()
    unmatched: List[str] = []
    for run in run_labels:
        hit = None
        for s in samples_cfg:
            sid = str(s.get("sample", ""))
            raw = s.get("raw_file")
            cand = {sid, _stem(sid), str(raw), _stem(raw)}
            if run in cand:
                hit = sid
                break
        if hit is None:  # casamento por substring como ultimo recurso
            for s in samples_cfg:
                sid = str(s.get("sample", ""))
                raw = str(s.get("raw_file", ""))
                if run and (run in sid or run in raw or _stem(raw) in run):
                    hit = sid
                    break
        if hit is None:
            unmatched.append(run)
        elif hit in used_samples:
            raise ValueError("sample %r casou com mais de uma coluna da matriz" % hit)
        else:
            rename[run] = hit
            used_samples.add(hit)
    if unmatched:
        raise ValueError(
            "Colunas da matriz sem mapeamento em samples.yaml: %s. "
            "Ajuste 'raw_file'/'sample' para casar com os nomes das colunas do motor." % unmatched
        )
    rows = [s for s in samples_cfg if str(s.get("sample")) in used_samples]
    df = pd.DataFrame(rows).set_index("sample")
    for c in ("condition", "bait_id", "role", "replicate", "raw_file"):
        if c not in df.columns:
            df[c] = None
    return rename, df[["condition", "bait_id", "role", "replicate", "raw_file"]]


def _bait_identifiers(baits_cfg: Optional[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Extrai, por isca, os identificadores possiveis (acessos/ids) para casar."""
    out: Dict[str, Dict[str, Any]] = {}
    if not baits_cfg:
        return out
    for b in baits_cfg.get("baits", []):
        bid = b.get("id")
        ids = set()
        if b.get("uniprot_accession"):
            ids.add(str(b["uniprot_accession"]))
        ids.add(str(bid))
        hdr = b.get("search_db_header") or ""
        # extrai token de acesso do cabecalho (ex.: 'bait|NSs_ICOARACI|...' ou 'sp|P17931|...')
        if "|" in hdr:
            parts = hdr.split("|")
            if len(parts) >= 2 and parts[1]:
                ids.add(parts[1])
        out[bid] = {
            "ids": {cont.strip_prefixes(x) for x in ids if x},
            "ortholog_accession": (b.get("ortholog") or {}).get("mouse_accession"),
            "source": b.get("source"),
        }
    return out


def _tag_baits(features: "pd.DataFrame", baits_cfg: Optional[Dict[str, Any]]) -> "pd.DataFrame":
    """Marca is_bait/species e sinaliza ortologo murino compartilhado com a isca."""
    info = _bait_identifiers(baits_cfg)
    if not info:
        return features
    all_ids = set()
    ortholog_accs = set()
    for meta in info.values():
        all_ids |= meta["ids"]
        if meta.get("ortholog_accession"):
            ortholog_accs.add(str(meta["ortholog_accession"]))

    for pid, row in features.iterrows():
        acc = str(row.get("leading_accession") or pid)
        majority = str(row.get("majority_protein_ids") or "")
        is_bait = acc in all_ids or any(x in majority for x in all_ids)
        if is_bait:
            features.at[pid, "is_bait"] = True
            features.at[pid, "species"] = "bait"
            # isca nunca e contaminante
            features.at[pid, "is_contaminant"] = False
        if acc in ortholog_accs:
            # ortologo murino da isca humana: prey legitimo, mas pode compartilhar peptideos
            features.at[pid, "shared_with_bait"] = True
    return features


def _assemble(
    intensity: "pd.DataFrame",
    features: "pd.DataFrame",
    run_labels: List[str],
    samples_cfg: Optional[List[Dict[str, Any]]],
    baits_cfg: Optional[Dict[str, Any]],
    provenance: Dict[str, Any],
) -> MSDataset:
    """Passos finais comuns: dedup por acesso, reconciliacao de amostras, iscas."""
    # dedup por acesso lider, mantendo a linha com mais peptides unicos
    if intensity.index.has_duplicates:
        order = features["n_unique_peptides"].fillna(0).astype(float).values
        pos = np.argsort(-order, kind="stable")   # maior evidencia primeiro
        intensity = intensity.iloc[pos]
        features = features.iloc[pos]
        keep = ~intensity.index.duplicated(keep="first")
        intensity = intensity[keep]
        features = features[keep]
        provenance = {**provenance, "deduped_accessions": True}

    rename, samples_df = _reconcile_samples(run_labels, samples_cfg)
    intensity = intensity.rename(columns=rename)
    # ordena colunas conforme samples_df
    intensity = intensity.reindex(columns=list(samples_df.index))
    features = _tag_baits(features, baits_cfg)
    ds = MSDataset(intensity=intensity, features=features, samples=samples_df, provenance=provenance)
    return ds.validate()


# ---------------------------------------------------------------------------
# MaxQuant (DDA)
# ---------------------------------------------------------------------------
def parse_maxquant(path, samples_cfg=None, baits_cfg=None) -> MSDataset:
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    # remove 'only identified by site' (pratica padrao Perseus)
    only_site = df.get("Only identified by site")
    n_only = int((only_site == "+").sum()) if only_site is not None else 0
    if only_site is not None:
        df = df[only_site != "+"].copy()

    id_col = "Majority protein IDs" if "Majority protein IDs" in df.columns else "Protein IDs"
    leading = df[id_col].map(_leading_accession)

    intens_cols, runs, rx = _match_intensity_columns(
        list(df.columns), [r"^LFQ intensity (.+)$", r"^Intensity (.+)$"]
    )
    if not intens_cols:
        raise ValueError("MaxQuant: nenhuma coluna 'LFQ intensity'/'Intensity' encontrada em %s" % path)
    quant_type = "lfq" if rx.startswith("^LFQ") else "intensity"

    intensity = pd.DataFrame(index=pd.Index(leading, name="protein_id"))
    for col, run in zip(intens_cols, runs):
        intensity[run] = _to_num(df[col].values, zero_to_nan=True)

    n = len(df)
    feats = empty_features(leading)
    feats["leading_accession"] = leading.values
    feats["majority_protein_ids"] = df[id_col].values
    feats["gene_names"] = _str_col(df, ["Gene names", "Gene Names"], n)
    feats["n_peptides"] = _num_col(df, "Peptides", n)
    feats["n_unique_peptides"] = _num_col(df, ["Razor + unique peptides", "Unique peptides"], n)
    if "Q-value" in df.columns:
        feats["q_value"] = pd.to_numeric(df["Q-value"], errors="coerce").values
    contam_col = "Potential contaminant" if "Potential contaminant" in df.columns else "Contaminant"
    contam_flag = (df.get(contam_col) == "+").values if contam_col in df.columns else np.zeros(n, bool)
    decoy_flag = (df.get("Reverse") == "+").values if "Reverse" in df.columns else np.zeros(n, bool)
    acc_contam, acc_decoy = cont.flag_series(leading.values)
    feats["is_contaminant"] = contam_flag | np.array(acc_contam)
    feats["is_decoy"] = decoy_flag | np.array(acc_decoy)

    prov = {
        "engine": "maxquant",
        "acquisition": "DDA",
        "quant_type": quant_type,
        "source_files": [str(path)],
        "n_only_identified_by_site_removed": n_only,
        "intensity_regex": rx,
    }
    return _assemble(intensity, feats, runs, samples_cfg, baits_cfg, prov)


# ---------------------------------------------------------------------------
# DIA-NN (DIA) - pg_matrix (+ report.tsv opcional)
# ---------------------------------------------------------------------------
def parse_diann(path, report_tsv=None, samples_cfg=None, baits_cfg=None) -> MSDataset:
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    id_cols = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"]
    present_id = [c for c in id_cols if c in df.columns]
    if "Protein.Group" not in df.columns:
        raise ValueError("DIA-NN: coluna 'Protein.Group' ausente em %s" % path)
    sample_cols = [c for c in df.columns if c not in id_cols]
    if not sample_cols:
        raise ValueError("DIA-NN: nenhuma coluna de amostra (MaxLFQ) encontrada em %s" % path)

    leading = df["Protein.Group"].map(_leading_accession)
    intensity = pd.DataFrame(index=pd.Index(leading, name="protein_id"))
    runs = []
    for col in sample_cols:
        run = Path(col).stem  # colunas costumam ser caminhos de arquivo
        runs.append(run)
        intensity[run] = _to_num(df[col].values, zero_to_nan=True)  # MaxLFQ ausente -> NaN

    feats = empty_features(leading)
    feats["leading_accession"] = leading.values
    feats["majority_protein_ids"] = df.get("Protein.Ids", df["Protein.Group"]).values
    feats["gene_names"] = df.get("Genes", pd.Series([""] * len(df))).values

    # evidencia de peptideos e q-value: do report.tsv se fornecido
    if report_tsv is not None and Path(report_tsv).exists():
        rep = pd.read_csv(report_tsv, sep="\t", dtype=str, low_memory=False)
        if {"Protein.Group", "Stripped.Sequence"}.issubset(rep.columns):
            npep = rep.groupby("Protein.Group")["Stripped.Sequence"].nunique()
            feats["n_unique_peptides"] = (
                pd.Series(feats["leading_accession"].values, index=feats.index)
                .map(lambda a: int(npep.get(a, 0)))
                .values
            )
            feats["n_peptides"] = feats["n_unique_peptides"].values
        if "PG.Q.Value" in rep.columns:
            qmin = pd.to_numeric(rep["PG.Q.Value"], errors="coerce").groupby(rep["Protein.Group"]).min()
            feats["q_value"] = (
                pd.Series(feats["leading_accession"].values, index=feats.index)
                .map(lambda a: qmin.get(a, np.nan))
                .values
            )

    acc_contam, acc_decoy = cont.flag_series(leading.values)
    feats["is_contaminant"] = np.array(acc_contam)
    feats["is_decoy"] = np.array(acc_decoy)

    prov = {
        "engine": "diann",
        "acquisition": "DIA",
        "quant_type": "maxlfq",
        "source_files": [str(path)] + ([str(report_tsv)] if report_tsv else []),
        "id_columns": present_id,
    }
    return _assemble(intensity, feats, runs, samples_cfg, baits_cfg, prov)


# ---------------------------------------------------------------------------
# FragPipe / Philosopher (DDA ou DIA)
# ---------------------------------------------------------------------------
def parse_fragpipe(path, samples_cfg=None, baits_cfg=None, acquisition="DDA") -> MSDataset:
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    acc_col = "Protein ID" if "Protein ID" in df.columns else "Protein"
    leading = df[acc_col].map(_leading_accession)

    intens_cols, runs, rx = _match_intensity_columns(
        list(df.columns), [r"(.+) MaxLFQ Intensity$", r"(.+) Intensity$"]
    )
    if not intens_cols:
        raise ValueError("FragPipe: nenhuma coluna '... MaxLFQ Intensity'/'... Intensity' em %s" % path)
    quant_type = "maxlfq" if "MaxLFQ" in rx else "intensity"
    zero_to_nan = acquisition.upper() == "DDA"

    intensity = pd.DataFrame(index=pd.Index(leading, name="protein_id"))
    for col, run in zip(intens_cols, runs):
        intensity[run] = _to_num(df[col].values, zero_to_nan=zero_to_nan)

    n = len(df)
    feats = empty_features(leading)
    feats["leading_accession"] = leading.values
    feats["majority_protein_ids"] = _str_col(df, ["Indistinguishable Proteins"], n) if "Indistinguishable Proteins" in df.columns else df[acc_col].values
    feats["gene_names"] = _str_col(df, ["Gene", "Gene Names"], n)
    feats["n_peptides"] = _num_col(df, ["Combined Total Peptides", "Total Peptides"], n)
    feats["n_unique_peptides"] = _num_col(df, ["Combined Unique Peptides", "Unique Peptides"], n)

    raw_protein = df.get("Protein", df[acc_col]).astype(str)
    acc_contam, acc_decoy = cont.flag_series(raw_protein.values)
    feats["is_contaminant"] = np.array(acc_contam)
    feats["is_decoy"] = np.array(acc_decoy)

    prov = {
        "engine": "fragpipe",
        "acquisition": acquisition.upper(),
        "quant_type": quant_type,
        "source_files": [str(path)],
        "intensity_regex": rx,
    }
    return _assemble(intensity, feats, runs, samples_cfg, baits_cfg, prov)


# ---------------------------------------------------------------------------
# Generico (CSV/TSV) - exige sidecar explicito
# ---------------------------------------------------------------------------
def parse_generic(path, sidecar: Dict[str, Any], samples_cfg=None, baits_cfg=None) -> MSDataset:
    """Matriz generica guiada por sidecar explicito (nunca adivinha).

    sidecar requer: id_col, e (intensity_cols OU intensity_regex). Opcionais:
    gene_col, sep (default '\t'), zero_to_nan (default True),
    acquisition (default 'DDA'), quant_type (default 'intensity').
    """
    path = Path(path)
    sep = sidecar.get("sep", "\t")
    df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)
    id_col = sidecar["id_col"]
    leading = df[id_col].map(_leading_accession)

    if sidecar.get("intensity_cols"):
        intens_cols = list(sidecar["intensity_cols"])
        runs = [str(c) for c in intens_cols]
    elif sidecar.get("intensity_regex"):
        intens_cols, runs, _ = _match_intensity_columns(list(df.columns), [sidecar["intensity_regex"]])
    else:
        raise ValueError("sidecar generico requer 'intensity_cols' ou 'intensity_regex'")
    if not intens_cols:
        raise ValueError("generico: nenhuma coluna de intensidade casou em %s" % path)

    zero_to_nan = bool(sidecar.get("zero_to_nan", True))
    intensity = pd.DataFrame(index=pd.Index(leading, name="protein_id"))
    for col, run in zip(intens_cols, runs):
        intensity[run] = _to_num(df[col].values, zero_to_nan=zero_to_nan)

    feats = empty_features(leading)
    feats["leading_accession"] = leading.values
    feats["majority_protein_ids"] = df[id_col].values
    if sidecar.get("gene_col") and sidecar["gene_col"] in df.columns:
        feats["gene_names"] = df[sidecar["gene_col"]].values
    acc_contam, acc_decoy = cont.flag_series(leading.values)
    feats["is_contaminant"] = np.array(acc_contam)
    feats["is_decoy"] = np.array(acc_decoy)

    prov = {
        "engine": "generic",
        "acquisition": str(sidecar.get("acquisition", "DDA")).upper(),
        "quant_type": sidecar.get("quant_type", "intensity"),
        "source_files": [str(path)],
    }
    return _assemble(intensity, feats, runs, samples_cfg, baits_cfg, prov)
