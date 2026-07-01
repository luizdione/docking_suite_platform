#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bonder.stats - estatistica quantitativa gold-standard (R) + QC (Python).

Inferencia (log2FC/p/FDR): DEP/limma/MSstats via runner.py (R/Bioconductor).
QC exploratorio e graficos: qc.py (Python), sem valor inferencial.
"""
from bonder.stats.prepare import write_dep_inputs
from bonder.stats.runner import run_dep, parse_dep_results
from bonder.stats import qc

__all__ = ["write_dep_inputs", "run_dep", "parse_dep_results", "qc"]
