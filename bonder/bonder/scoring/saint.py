#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""saint.py - integracao com SAINTexpress (Teo et al., J Proteomics 2014).

SAINTexpress e o padrao da comunidade para escore de especificidade isca-presa
em AP-MS/CoIP-MS: da a probabilidade SAINT e o BFDR de cada prey ser interator
verdadeiro, comparando os IPs da isca contra os controles.

NAO reimplementamos o algoritmo. Aqui:
  - formatamos os 3 arquivos de entrada exigidos (interaction/prey/bait);
  - chamamos o binario SAINTexpress-int (versao para intensidades, ex.: LFQ/MaxLFQ);
  - parseamos a saida 'list.txt' num DataFrame.

O binario SAINTexpress roda dentro do conteiner de scoring (Docker); em CPU sem o
binario, write_saint_inputs continua funcionando (util para testes e para gerar
os arquivos para rodar em outra maquina).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["numpy", "pandas"])
    import numpy as np
    import pandas as pd

DEFAULT_PREY_LENGTH = 500  # usado quando o comprimento real da proteina e desconhecido


def write_saint_inputs(
    ds,
    contrast: Dict[str, Any],
    out_dir,
    prey_lengths: Optional[Dict[str, int]] = None,
    gene_map: Optional[Dict[str, str]] = None,
    exclude_baits: bool = True,
) -> Dict[str, Path]:
    """Gera interaction.txt, prey.txt e bait.txt para SAINTexpress-int.

    contrast: dict com 'name', 'bait_id', 'test' (condicao IP) e 'control'.
    Retorna dict com os caminhos dos 3 arquivos.

    Formatos (tab-separado, sem cabecalho):
      interaction : IP <TAB> Bait <TAB> Prey <TAB> Intensity
      prey        : Prey <TAB> PreyLength <TAB> PreyGene
      bait        : IP <TAB> Bait <TAB> {T|C}
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prey_lengths = prey_lengths or {}
    gene_map = gene_map or {}

    bait_id = contrast["bait_id"]
    test_cond = contrast["test"]
    ctrl_cond = contrast["control"]
    test_ips = ds.samples_in(test_cond)
    ctrl_ips = ds.samples_in(ctrl_cond)
    if not test_ips or not ctrl_ips:
        raise ValueError(
            "Contraste %r sem replicas: test=%d control=%d"
            % (contrast.get("name"), len(test_ips), len(ctrl_ips))
        )

    # preys candidatos = todas as proteinas, opcionalmente sem as iscas exogenas
    feats = ds.features
    prey_index = list(ds.intensity.index)
    if exclude_baits:
        bait_rows = set(feats.index[feats["is_bait"].fillna(False).astype(bool)])
        prey_index = [p for p in prey_index if p not in bait_rows]

    # --- interaction file ---
    inter_rows: List[str] = []
    for ip in test_ips + ctrl_ips:
        col = ds.intensity[ip]
        for prey in prey_index:
            val = col.get(prey, np.nan)
            if pd.notna(val) and val > 0:
                inter_rows.append("%s\t%s\t%s\t%g" % (ip, bait_id, prey, float(val)))
    inter_path = out_dir / "interaction.txt"
    inter_path.write_text("\n".join(inter_rows) + "\n", encoding="utf-8")

    # --- prey file ---
    prey_rows: List[str] = []
    for prey in prey_index:
        length = int(prey_lengths.get(prey, DEFAULT_PREY_LENGTH))
        gene = gene_map.get(prey) or str(feats.at[prey, "gene_names"] or prey)
        gene = gene.split(";")[0] if gene else prey
        prey_rows.append("%s\t%d\t%s" % (prey, length, gene))
    prey_path = out_dir / "prey.txt"
    prey_path.write_text("\n".join(prey_rows) + "\n", encoding="utf-8")

    # --- bait file --- (T = test/isca ; C = controle)
    bait_rows_out = [("%s\t%s\tT" % (ip, bait_id)) for ip in test_ips]
    bait_rows_out += [("%s\t%s\tC" % (ip, bait_id)) for ip in ctrl_ips]
    bait_path = out_dir / "bait.txt"
    bait_path.write_text("\n".join(bait_rows_out) + "\n", encoding="utf-8")

    return {"interaction": inter_path, "prey": prey_path, "bait": bait_path}


def run_saintexpress(
    inputs: Dict[str, Path],
    out_dir,
    binary: str = "SAINTexpress-int",
    extra_args: Optional[List[str]] = None,
) -> Path:
    """Executa SAINTexpress-int nos arquivos gerados. Retorna o caminho de list.txt.

    Levanta FileNotFoundError com instrucao clara se o binario nao existir
    (ele mora no conteiner de scoring; ver bonder/docker/Dockerfile.scoring).
    """
    out_dir = Path(out_dir)
    exe = shutil.which(binary)
    if exe is None:
        raise FileNotFoundError(
            "Binario %r nao encontrado. SAINTexpress roda no conteiner de scoring; "
            "para rodar localmente compile de https://github.com/rachelpetit/SAINTexpress "
            "ou use 'docker compose run scoring'. Os inputs ja foram gerados em %s."
            % (binary, out_dir)
        )
    cmd = [exe, str(inputs["interaction"]), str(inputs["prey"]), str(inputs["bait"])]
    if extra_args:
        cmd += list(extra_args)
    subprocess.run(cmd, cwd=str(out_dir), check=True)
    list_path = out_dir / "list.txt"
    if not list_path.exists():
        raise RuntimeError("SAINTexpress terminou sem gerar list.txt em %s" % out_dir)
    return list_path


def parse_saint_output(list_path) -> "pd.DataFrame":
    """Le o 'list.txt' do SAINTexpress num DataFrame padronizado.

    Colunas de interesse: Bait, Prey, PreyGene, AvgSpec (AvgInt), FoldChange,
    SaintScore, BFDR. Renomeadas para minusculas estaveis.
    """
    df = pd.read_csv(list_path, sep="\t", dtype=str)
    rename = {
        "Bait": "bait", "Prey": "prey", "PreyGene": "prey_gene",
        "AvgSpec": "avg_intensity", "AvgInt": "avg_intensity",
        "FoldChange": "fold_change", "SaintScore": "saint_score", "BFDR": "bfdr",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    for c in ("avg_intensity", "fold_change", "saint_score", "bfdr"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df
