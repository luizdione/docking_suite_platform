#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
md_clean.py  —  SHIM DE COMPATIBILIDADE (DEPRECIADO)

O pipeline de Dinamica Molecular foi RENOMEADO para `md_simulation.py`
(refatoracao de 2026-05-25). Este arquivo continua funcionando apenas para
nao quebrar comandos antigos, .bat, scripts e referencias historicas:
ele delega 100% para md_simulation.main().

>>> Atualize seus comandos para:
        python md_simulation.py --analyze-only --tag prod_v3
        python md_simulation.py --fresh --ns 50 --tag prod_v4 --platform CUDA ...

(O diretorio de saida continua sendo docking_project/figures/md_clean/ —
 so o SCRIPT mudou de nome, nao a pasta de resultados.)
"""
import sys

print("[aviso] 'md_clean.py' foi renomeado para 'md_simulation.py'. "
      "Delegando a execucao; atualize seus comandos.", file=sys.stderr)

from md_simulation import main  # noqa: E402  (mantem a API: tudo vive em md_simulation)

if __name__ == "__main__":
    main()
