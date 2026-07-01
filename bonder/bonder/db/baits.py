#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""baits.py - registro de iscas + montagem do FASTA de busca combinado.

As iscas do piloto sao EXOGENAS ao proteoma murino. Para que o motor de busca
consiga quantificar cada isca (controle positivo), o FASTA de busca deve conter:
  proteoma murino (UniProt 10090)  +  contaminantes (cRAP)  +  sequencias das iscas.

Este modulo le config/baits.yaml, valida as sequencias das iscas e concatena o
search_db.fasta com cabecalhos estaveis (search_db_header). A isca humana LGALS3
(P17931) e mantida DISTINTA do ortologo murino Lgals3 (P16110).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from bonder.db.uniprot import iter_fasta


def _read_single_fasta_seq(path: Path) -> str:
    seqs = [seq for _acc, _h, seq in iter_fasta(path)]
    if not seqs:
        raise ValueError("FASTA de isca vazio ou invalido: %s" % path)
    return seqs[0]


def resolve_bait_sequences(baits_cfg: Dict[str, Any], config_dir: Path) -> Dict[str, str]:
    """id da isca -> sequencia de aminoacidos (lida do FASTA apontado em baits.yaml)."""
    out: Dict[str, str] = {}
    for b in baits_cfg.get("baits", []):
        fasta_rel = b.get("fasta")
        if not fasta_rel:
            continue
        fasta_path = (config_dir / fasta_rel).resolve()
        if fasta_path.exists():
            out[b["id"]] = _read_single_fasta_seq(fasta_path)
    return out


def build_search_fasta(
    baits_cfg: Dict[str, Any],
    config_dir: Path,
    out_path: Path,
    mouse_fasta: Optional[Path] = None,
    contaminants_fasta: Optional[Path] = None,
) -> Path:
    """Concatena proteoma murino + contaminantes + iscas em um search_db.fasta.

    mouse_fasta / contaminants_fasta sao opcionais (podem faltar em ambiente de
    teste); nesse caso o FASTA gerado contem ao menos as iscas, permitindo que o
    usuario complete com o proteoma baixado por 'bonder db sync'.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    seqs = resolve_bait_sequences(baits_cfg, config_dir)

    def _wrap(seq: str, width: int = 60) -> str:
        return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

    with open(out_path, "w", encoding="utf-8") as out:
        for src in (mouse_fasta, contaminants_fasta):
            if src and Path(src).exists():
                for _acc, header, seq in iter_fasta(src):
                    out.write(">%s\n%s\n" % (header, _wrap(seq)))
        for b in baits_cfg.get("baits", []):
            bid = b["id"]
            if bid in seqs:
                header = b.get("search_db_header") or ("bait|%s|%s" % (bid, b.get("display_name", bid)))
                out.write(">%s\n%s\n" % (header, _wrap(seqs[bid])))
    return out_path


def validate_baits(baits_cfg: Dict[str, Any], config_dir: Path) -> List[str]:
    """Retorna lista de avisos (iscas sem FASTA, ortologo nao declarado, etc.)."""
    warnings: List[str] = []
    for b in baits_cfg.get("baits", []):
        bid = b.get("id")
        fasta_rel = b.get("fasta")
        if not fasta_rel or not (config_dir / fasta_rel).exists():
            warnings.append(
                "Isca %s sem FASTA em %s: forneca a sequencia para o controle positivo "
                "e o FASTA de busca." % (bid, fasta_rel)
            )
        if b.get("source") == "human" and not (b.get("ortholog") or {}).get("mouse_accession"):
            warnings.append("Isca humana %s sem ortologo murino declarado (ortholog.mouse_accession)." % bid)
    return warnings
