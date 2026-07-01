#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""prepare.py - gera entradas para a estatistica gold-standard em R.

Para matrizes em nivel de PROTEINA (proteinGroups/pg_matrix/combined_protein) a
ferramenta de referencia e o DEP (Zhang et al., Nat Protoc 2018), que envolve o
limma com filtragem, normalizacao (vsn) e imputacao MNAR (MinProb) apropriadas
para enriquecimento AP-MS. Quando ha dados em nivel de FEATURE, o MSstats e
preferivel; aqui exportamos o formato que o script R correspondente consome.

NAO reimplementamos o teste: apenas escrevemos os arquivos.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["pandas"])
    import pandas as pd


def write_dep_inputs(ds, contrast: Dict[str, Any], out_dir) -> Dict[str, Path]:
    """Escreve data.tsv (wide) e experimental_design.tsv para run_dep_limma.R.

    data.tsv: colunas ID, name (gene), e uma coluna de intensidade por amostra.
    experimental_design.tsv: label (=coluna de intensidade), condition, replicate.
    Restringe as amostras do contraste (test + control).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    test_cond, ctrl_cond = contrast["test"], contrast["control"]
    samples = ds.samples_in(test_cond) + ds.samples_in(ctrl_cond)
    if not samples:
        raise ValueError("Contraste %r sem amostras." % contrast.get("name"))

    intens = ds.intensity[samples].copy()
    genes = ds.features["gene_names"].reindex(intens.index).fillna("")
    data = intens.copy()
    data.insert(0, "name", [g.split(";")[0] if isinstance(g, str) and g else acc
                            for acc, g in zip(intens.index, genes)])
    data.insert(0, "ID", intens.index)
    data_path = out_dir / "data.tsv"
    data.to_csv(data_path, sep="\t", index=False)

    smp = ds.samples.loc[samples]
    design = pd.DataFrame({
        "label": samples,
        "condition": smp["condition"].values,
        "replicate": smp["replicate"].values,
    })
    design_path = out_dir / "experimental_design.tsv"
    design.to_csv(design_path, sep="\t", index=False)

    return {"data": data_path, "design": design_path,
            "test": test_cond, "control": ctrl_cond}
