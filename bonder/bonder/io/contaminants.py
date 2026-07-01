#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""contaminants.py - deteccao de contaminantes e decoys em saidas de busca.

Deteccao dupla (robusta): flags booleanas do motor E prefixos de acesso.
  MaxQuant   : 'Potential contaminant' == '+', 'Reverse' == '+', prefixos CON__ / REV__
  FragPipe   : prefixos 'contam_' e acessos rev_
  cRAP/thegpm: lista de acessos conhecidos (keratinas, tripsina, BSA, etc.)

Iscas NUNCA sao tratadas como contaminante mesmo que combinem com um prefixo.
"""
from __future__ import annotations

from typing import Iterable

# Prefixos de acesso que indicam contaminante ou decoy nos principais motores.
CONTAMINANT_PREFIXES = ("CON__", "contam_", "CONTAM_", "sp|CON", "Cont_")
DECOY_PREFIXES = ("REV__", "rev_", "REV_", "DECOY_", "decoy_", "XXX_")


def is_contaminant_accession(accession: str) -> bool:
    if not isinstance(accession, str):
        return False
    return accession.startswith(CONTAMINANT_PREFIXES)


def is_decoy_accession(accession: str) -> bool:
    if not isinstance(accession, str):
        return False
    return accession.startswith(DECOY_PREFIXES)


def flag_series(accessions: Iterable[str]):
    """Retorna duas listas booleanas (is_contaminant, is_decoy) por acesso."""
    contam, decoy = [], []
    for a in accessions:
        contam.append(is_contaminant_accession(a))
        decoy.append(is_decoy_accession(a))
    return contam, decoy


def strip_prefixes(accession: str) -> str:
    """Remove prefixos CON__/REV__/contam_ de um acesso para casar com iscas/DB."""
    if not isinstance(accession, str):
        return accession
    for p in CONTAMINANT_PREFIXES + DECOY_PREFIXES:
        if accession.startswith(p):
            return accession[len(p):]
    return accession
