#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Testes de formatacao de entrada do SAINTexpress (CPU, sem o binario)."""
import pandas as pd

from bonder.io.parsers import parse_maxquant
from bonder.scoring import write_saint_inputs
from tests.test_io_parsers import _write_maxquant, SAMPLES_CFG, BAITS_CFG

CONTRAST = {"name": "NSs_vs_CTRL", "bait_id": "NSs_ICOARACI", "test": "IP_NSs", "control": "CTRL_NSs"}


def test_saint_inputs_shape(tmp_path):
    p = tmp_path / "proteinGroups.txt"
    _write_maxquant(p)
    ds = parse_maxquant(p, samples_cfg=SAMPLES_CFG, baits_cfg=BAITS_CFG)
    outdir = tmp_path / "saint"
    paths = write_saint_inputs(ds, CONTRAST, outdir, exclude_baits=True)

    bait = paths["bait"].read_text().strip().splitlines()
    # 2 IPs teste (T) + 2 controles (C)
    assert sum(l.endswith("\tT") for l in bait) == 2
    assert sum(l.endswith("\tC") for l in bait) == 2

    inter = paths["interaction"].read_text().strip().splitlines()
    # cada linha: IP, Bait, Prey, Intensity
    assert all(len(l.split("\t")) == 4 for l in inter)
    # a isca exogena NAO deve aparecer como prey (exclude_baits=True)
    assert all(l.split("\t")[2] != "NSs_ICOARACI" for l in inter)
    # apenas intensidades > 0 entram
    assert all(float(l.split("\t")[3]) > 0 for l in inter)

    prey = paths["prey"].read_text().strip().splitlines()
    assert all(len(l.split("\t")) == 3 for l in prey)


def test_saint_parse_output(tmp_path):
    from bonder.scoring import parse_saint_output
    lst = tmp_path / "list.txt"
    pd.DataFrame({
        "Bait": ["NSs_ICOARACI", "NSs_ICOARACI"],
        "Prey": ["P12345", "Q8R1"],
        "PreyGene": ["Prey1", "Prey2"],
        "AvgInt": ["1050", "500"],
        "FoldChange": ["10.5", "2.1"],
        "SaintScore": ["0.99", "0.40"],
        "BFDR": ["0.00", "0.15"],
    }).to_csv(lst, sep="\t", index=False)
    out = parse_saint_output(lst)
    assert "saint_score" in out.columns and "bfdr" in out.columns
    assert out.loc[0, "saint_score"] == 0.99
