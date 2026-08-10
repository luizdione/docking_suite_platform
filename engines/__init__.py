# -*- coding: utf-8 -*-
"""
Pacote engines/ -- abstracao de motor da plataforma (Fase 4).

Motores: vina (docking), boltz2 (estrutura+afinidade), openmm (MD), haddock3
(pose por restricoes). Todos implementam o contrato de engines.base.Engine e se
registram via @register; o runner.py resolve o engine de cada task por nome.

Imports pesados (docking_v2, md_simulation, openmm, boltz) sao TARDIOS, dentro
de check_toolchain()/run() de cada engine -- importar `engines` continua leve.
"""
from engines.base import (
    Engine,
    EngineRequest,
    EngineResult,
    get_engine,
    list_engines,
    register,
)
from engines import toolchain

__all__ = [
    "Engine",
    "EngineRequest",
    "EngineResult",
    "get_engine",
    "list_engines",
    "register",
    "toolchain",
]
