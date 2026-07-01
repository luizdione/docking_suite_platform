#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mist.py - formatacao de entrada para MiST (Verschueren et al. / Jager et al.).

MiST (Mass spectrometry interaction STatistics) combina tres metricas por prey -
reprodutibilidade, abundancia e especificidade - num escore unico. E o padrao
usado em interactomas host-pathogen (ex.: interactoma SARS-CoV-2, Gordon et al.,
Nature 2020), por isso e a segunda pontuacao recomendada para a isca VIRAL NSs.

Aqui geramos a matriz de entrada no layout esperado pelo MiST (preys x experimentos
com um mapa de chaves bait/experimento). A execucao do MiST em si e feita pela
implementacao publicada dentro do conteiner (via artMS/MiST no R, ou o script
MiST original apontado por BONDER_MIST_CMD). NAO reimplementamos o escore.
"""
from __future__ import annotations

import os
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


def write_mist_inputs(ds, contrast: Dict[str, Any], out_dir, exclude_baits: bool = True) -> Dict[str, Path]:
    """Gera data matrix (preys x IPs) e keys (IP -> bait/condicao) para MiST.

    data.txt : primeira coluna 'Prey' + uma coluna por IP (intensidade; 0 se ausente)
    keys.txt : Experiment <TAB> Bait <TAB> Condition (bait para IPs, control p/ ctrls)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bait_id = contrast["bait_id"]
    test_ips = ds.samples_in(contrast["test"])
    ctrl_ips = ds.samples_in(contrast["control"])
    ips = test_ips + ctrl_ips

    feats = ds.features
    idx = list(ds.intensity.index)
    if exclude_baits:
        bait_rows = set(feats.index[feats["is_bait"].fillna(False).astype(bool)])
        idx = [p for p in idx if p not in bait_rows]

    data = ds.intensity.loc[idx, ips].fillna(0.0)
    data.index.name = "Prey"
    data_path = out_dir / "data.txt"
    data.reset_index().to_csv(data_path, sep="\t", index=False)

    keys = pd.DataFrame({
        "Experiment": ips,
        "Bait": [bait_id] * len(ips),
        "Condition": (["bait"] * len(test_ips)) + (["control"] * len(ctrl_ips)),
    })
    keys_path = out_dir / "keys.txt"
    keys.to_csv(keys_path, sep="\t", index=False)

    return {"data": data_path, "keys": keys_path}


def run_mist(inputs: Dict[str, Path], out_dir) -> Optional[Path]:
    """Executa a implementacao publicada do MiST se BONDER_MIST_CMD estiver definido.

    BONDER_MIST_CMD deve ser um comando que aceite: <data> <keys> <out>. Sem ele,
    retorna None (o pipeline segue com SAINT/MSstats e registra a ausencia).
    """
    cmd_tmpl = os.environ.get("BONDER_MIST_CMD")
    if not cmd_tmpl:
        return None
    out_dir = Path(out_dir)
    out_path = out_dir / "mist_scores.txt"
    parts = cmd_tmpl.split() + [str(inputs["data"]), str(inputs["keys"]), str(out_path)]
    subprocess.run(parts, cwd=str(out_dir), check=True)
    return out_path if out_path.exists() else None
