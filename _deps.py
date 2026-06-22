#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_deps.py — politica centralizada de dependencias da plataforma NQ-Docking Suite.

Reprodutibilidade primeiro: por padrao a plataforma NAO instala pacotes em
runtime (anti-padrao para ambientes controlados/CI — mascara um env quebrado e
torna o resultado nao-reproduzivel). Quando uma dependencia falta, a funcao
`ensure()` levanta ImportError com uma instrucao clara apontando o requirements.

A auto-instalacao continua disponivel, mas apenas como OPT-IN explicito:
    Windows : set DOCKING_AUTO_PIP=1
    bash    : export DOCKING_AUTO_PIP=1

Uso tipico (no bloco except ImportError de um modulo):
    try:
        import numpy as np
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        import _deps
        _deps.ensure(["numpy", "matplotlib"])  # raise, ou instala se opt-in
        import numpy as np
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
"""
import os
import sys
import subprocess

#: True somente se o usuario habilitou explicitamente a auto-instalacao.
AUTO_PIP = os.environ.get("DOCKING_AUTO_PIP", "").strip().lower() in ("1", "true", "yes", "on")


def ensure(packages, *, hint=None):
    """Garante dependencias ausentes — sem instalar nada silenciosamente.

    Chamar dentro de um `except ImportError`. Comportamento:
      - DOCKING_AUTO_PIP habilitado  -> instala `packages` via pip e retorna True
        (o chamador reexecuta os imports na sequencia).
      - caso contrario               -> levanta ImportError com mensagem clara
        apontando requirements.txt (reprodutibilidade).

    packages: str ou iteravel de specs pip (ex.: ["numpy", "matplotlib"]).
    hint: texto extra opcional para a mensagem de erro.
    """
    if isinstance(packages, str):
        packages = [packages]
    packages = [str(p) for p in packages]
    if AUTO_PIP:
        subprocess.run([sys.executable, "-m", "pip", "install", *packages, "--quiet"],
                       check=False)
        return True
    msg = ("Dependencia(s) ausente(s): " + ", ".join(packages) +
           ". Instale o ambiente do projeto: 'pip install -r requirements.txt' "
           "(ou use o conda env 'datascience'). "
           "Para auto-instalar em runtime defina a variavel DOCKING_AUTO_PIP=1.")
    if hint:
        msg += " " + hint
    raise ImportError(msg)
