# -*- coding: utf-8 -*-
"""
engines/base.py -- contrato comum dos motores (Fase 4).

Um Engine e a fronteira entre "o que rodar" (RunSpec/plan do runner, que ja
existe) e "como rodar" (Vina, Boltz-2, OpenMM, HADDOCK3). O runner.py da Fase 6
foi escrito com um executor INJETAVEL exatamente para este momento: aqui a
implementacao default de cada engine entra sem reescrever o runner.

Contrato (Engine ABC):
    capabilities()    -> set[str]     o que o engine entrega: {"pose","affinity",
                                       "structure","confidence","md"}
    check_toolchain() -> dict         {available: bool, tools: {...}} -- resolve
                                       binarios SEM instalar nada (ver toolchain.py)
    prepare(req)      -> EngineRequest resolve entradas/saidas; nao executa
    run(req)          -> EngineResult  executa (ou planeja, se req.dry_run)
    parse(result,req) -> dict          metricas normalizadas

Regra de ouro herdada do resto da plataforma: NUNCA degradar em silencio. Falta
de toolchain ou de entrada vira EngineResult(status="error", error=...) com
mensagem clara -- nunca um resultado plausivel e errado. Gates de QC (torsoes,
RMSD, pLDDT) NAO moram aqui: sao aplicados pelo runner sobre os artefatos, para
que "reprovar num gate" continue sendo estado de dados, nao excecao.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class EngineRequest:
    """Tudo o que um engine precisa para uma etapa de um par alvo-ligante."""
    campaign: str
    target: str
    ligand: str
    kind: str                              # docking|structure|md|pose
    protocol: dict = field(default_factory=dict)   # protocolo resolvido do registry
    site: str | None = None
    params: dict = field(default_factory=dict)     # grid_center, exhaustiveness, seed...
    inputs: dict = field(default_factory=dict)     # role -> path (receptor_pdbqt, ...)
    outputs: dict = field(default_factory=dict)    # role -> path (preenchido em prepare/run)
    target_info: dict = field(default_factory=dict)  # projecao do registry (pdb, sites...)
    workdir: Path | None = None
    dry_run: bool = False


@dataclass
class EngineResult:
    """Resultado da EXECUCAO (nao do QC). status: ok|error|skipped."""
    engine: str
    status: str = "ok"
    cmdline: list | None = None
    outputs: dict = field(default_factory=dict)
    metrics: dict = field(default_factory=dict)
    log: str = ""
    error: str | None = None

    @property
    def ok(self):
        return self.status == "ok" and not self.error


class Engine(ABC):
    """Contrato de um motor. Subclasses declaram `name` (str)."""

    name: str = "base"

    @abstractmethod
    def capabilities(self) -> set:
        """Conjunto de capacidades: {'pose','affinity','structure','confidence','md'}."""

    @abstractmethod
    def check_toolchain(self) -> dict:
        """{available: bool, tools: {nome: {...}}}. Nao instala nada; so resolve."""

    def prepare(self, req: EngineRequest) -> EngineRequest:
        """Resolve entradas/saidas. Default: no-op (engine simples). Nao executa."""
        return req

    @abstractmethod
    def run(self, req: EngineRequest) -> EngineResult:
        """Executa (ou planeja, se req.dry_run). NAO levanta por erro de dominio:
        devolve EngineResult(status='error', error=...)."""

    def parse(self, result: EngineResult, req: EngineRequest) -> dict:
        """Metricas normalizadas do resultado. Default: as proprias metrics."""
        return dict(result.metrics)

    def __repr__(self):
        return f"<Engine {self.name}>"


# --------------------------------------------------------------- registro

_ENGINES = {}


def register(cls):
    """Decorator: registra uma subclasse de Engine pelo seu `name`."""
    name = getattr(cls, "name", None)
    if not name or name == "base":
        raise ValueError(f"engine {cls!r} precisa de um `name` proprio")
    _ENGINES[name] = cls
    return cls


def get_engine(name: str) -> Engine:
    """Instancia o engine registrado com este nome. Erra alto se nao existe."""
    _ensure_loaded()
    try:
        return _ENGINES[name]()
    except KeyError:
        raise KeyError(
            f"engine '{name}' nao registrado. Disponiveis: {sorted(_ENGINES)}"
        ) from None


def list_engines() -> list:
    _ensure_loaded()
    return sorted(_ENGINES)


_LOADED = False


def _ensure_loaded():
    """Importa os modulos de engine sob demanda, para registra-los. Import tardio
    e tolerante: um engine com dependencia pesada quebrada nao derruba os demais
    (a indisponibilidade aparece depois em check_toolchain)."""
    global _LOADED
    if _LOADED:
        return
    _LOADED = True
    from importlib import import_module
    for mod in ("vina", "boltz2", "openmm_md", "haddock3"):
        try:
            import_module(f"engines.{mod}")
        except Exception:  # noqa: BLE001 -- registro tardio nao pode ser fatal
            pass
