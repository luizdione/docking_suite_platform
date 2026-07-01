#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""esm2.py - priorizacao EXPLORATORIA por embeddings ESM2 (SUPLEMENTAR).

ATENCAO cientifica: embeddings ESM2 NAO sao metodo padrao para CHAMAR interatores.
Este modulo entra apenas como camada exploratoria de re-ranqueamento/priorizacao,
NUNCA como base do escore de confianca que vai para a figura principal (esse vem
de SAINT/DEP). Fica desligado por padrao (config/gpu.yaml).

Boas praticas implementadas: chunk <= max_residues, mean-pool com mascara de
atencao em fp32, L2-normalizacao, cache em disco por (sha256 da seq, modelo). Sem
torch/transformers, degrada com mensagem clara.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional

from bonder.gpu.detect import probe


def _seq_hash(seq: str, model: str) -> str:
    return hashlib.sha256(("%s|%s" % (model, seq)).encode("utf-8")).hexdigest()


def embed_sequences(
    sequences: Dict[str, str],
    model_gpu: str = "facebook/esm2_t33_650M_UR50D",
    model_cpu: str = "facebook/esm2_t30_150M_UR50D",
    max_residues: int = 1022,
    batch_size: int = 8,
    cache_dir: Optional[Path] = None,
):
    """accession -> vetor de embedding (numpy float32, L2-normalizado).

    Retorna dict {accession: np.ndarray}. Levanta ImportError se torch/transformers
    faltarem (o chamador trata como 'priorizacao ESM2 indisponivel').
    """
    try:
        import numpy as np
        import torch
        from transformers import AutoTokenizer, AutoModel
    except ImportError as e:
        raise ImportError(
            "ESM2 requer torch + transformers (conteiner GPU). Priorizacao exploratoria "
            "indisponivel: %s" % e
        )

    dec = probe()
    device = dec.device
    model_name = model_gpu if device == "cuda" else model_cpu
    dtype = torch.float16 if device == "cuda" else torch.float32

    cache: Dict[str, "np.ndarray"] = {}
    if cache_dir:
        cache_dir = Path(cache_dir); cache_dir.mkdir(parents=True, exist_ok=True)

    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    if device == "cuda":
        model = model.half()
    model = model.to(device).eval()

    out: Dict[str, "np.ndarray"] = {}
    items = list(sequences.items())
    for i in range(0, len(items), batch_size):
        batch = items[i:i + batch_size]
        ids, seqs = [], []
        for acc, seq in batch:
            if cache_dir:
                cf = cache_dir / ("%s.npy" % _seq_hash(seq, model_name))
                if cf.exists():
                    out[acc] = np.load(cf)
                    continue
            ids.append(acc)
            seqs.append(seq[:max_residues])
        if not seqs:
            continue
        enc = tok(seqs, return_tensors="pt", padding=True, truncation=True, max_length=max_residues)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            rep = model(**enc).last_hidden_state  # (B, L, H)
        mask = enc["attention_mask"].unsqueeze(-1).to(torch.float32)
        summed = (rep.to(torch.float32) * mask).sum(dim=1)
        counts = mask.sum(dim=1).clamp(min=1.0)
        pooled = summed / counts
        norm = pooled / pooled.norm(dim=1, keepdim=True).clamp(min=1e-8)
        vecs = norm.cpu().numpy().astype("float32")
        for acc, v in zip(ids, vecs):
            out[acc] = v
            if cache_dir:
                np.save(cache_dir / ("%s.npy" % _seq_hash(sequences[acc], model_name)), v)
        if device == "cuda":
            torch.cuda.empty_cache()
    return out


def cosine_prioritize(embeddings, reference_accessions: List[str]):
    """Escore de plausibilidade = cosseno medio vs conjunto de interatores conhecidos."""
    import numpy as np
    refs = [embeddings[a] for a in reference_accessions if a in embeddings]
    if not refs:
        return {}
    ref = np.mean(np.stack(refs), axis=0)
    ref = ref / (np.linalg.norm(ref) + 1e-8)
    return {acc: float(np.dot(v, ref)) for acc, v in embeddings.items()}
