#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bonder.scoring - especificidade isca-presa via ferramentas consagradas.

SAINTexpress (primario, todas as iscas) e MiST (host-pathogen, isca viral NSs).
Nenhum algoritmo e reimplementado: formatamos entradas e envolvemos os binarios
publicados.
"""
from bonder.scoring.saint import (
    write_saint_inputs,
    run_saintexpress,
    parse_saint_output,
)
from bonder.scoring.mist import write_mist_inputs, run_mist

__all__ = [
    "write_saint_inputs",
    "run_saintexpress",
    "parse_saint_output",
    "write_mist_inputs",
    "run_mist",
]
