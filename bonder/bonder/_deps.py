#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""_deps.py - politica centralizada de dependencias do bonder.

Reprodutibilidade primeiro: por padrao a plataforma NAO instala pacotes em
runtime (anti-padrao para ambientes controlados/CI - mascara um env quebrado e
torna o resultado nao-reproduzivel). Quando uma dependencia falta, `ensure()`
levanta ImportError com instrucao clara apontando o requirements/ambiente conda.

Auto-instalacao existe apenas como OPT-IN explicito (mesma convencao da suite):
    Windows : set BONDER_AUTO_PIP=1
    bash    : export BONDER_AUTO_PIP=1

Uso tipico (no bloco except ImportError de um modulo):
    try:
        import pandas as pd
    except ImportError:
        from bonder import _deps
        _deps.ensure(["pandas"])
        import pandas as pd
"""
import os
import sys
import subprocess

#: True somente se o usuario habilitou explicitamente a auto-instalacao.
AUTO_PIP = os.environ.get("BONDER_AUTO_PIP", "").strip().lower() in ("1", "true", "yes", "on")


def ensure(packages, *, hint=None):
    """Garante dependencias ausentes - sem instalar nada silenciosamente.

    Chamar dentro de um `except ImportError`. Comportamento:
      - BONDER_AUTO_PIP habilitado -> instala `packages` via pip e retorna True
        (o chamador reexecuta os imports na sequencia).
      - caso contrario -> levanta ImportError com mensagem clara.

    packages: str ou iteravel de specs pip (ex.: ["pandas", "pyyaml"]).
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
           ". Instale o ambiente do projeto: 'pip install -r bonder/requirements.txt' "
           "ou use o conda env 'datascience' (conda env update -f bonder/environment.yml). "
           "Para auto-instalar em runtime defina BONDER_AUTO_PIP=1.")
    if hint:
        msg += " " + hint
    raise ImportError(msg)
