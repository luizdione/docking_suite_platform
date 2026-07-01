#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""config.py - carregamento centralizado de configuracao do bonder.

Le os YAMLs em bonder/config/ (baits, samples, databases, gpu), resolve caminhos
e aplica sobrescritas por variavel de ambiente:

    BONDER_CONFIG_DIR   diretorio dos YAMLs (default: <repo>/config)
    BONDER_DATA_ROOT    raiz do cache de dados de referencia
    BONDER_OFFLINE      '1' -> modo offline (usa so o cache)
    BONDER_RUNS_DIR     diretorio de saida das corridas (default: <repo>/runs)

Tudo e ASCII-safe e nao instala nada em runtime.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import yaml
except ImportError:  # pragma: no cover - politica de dependencia
    from bonder import _deps
    _deps.ensure(["pyyaml"])
    import yaml


# Raiz do subprojeto bonder = pai do pacote 'bonder'.
PACKAGE_DIR = Path(__file__).resolve().parent          # .../bonder/bonder
PROJECT_ROOT = PACKAGE_DIR.parent                       # .../bonder


def _config_dir() -> Path:
    env = os.environ.get("BONDER_CONFIG_DIR")
    return Path(env).expanduser().resolve() if env else (PROJECT_ROOT / "config")


def _load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(
            "Config ausente: %s. Verifique BONDER_CONFIG_DIR ou o diretorio config/." % path
        )
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError("Config invalida (esperado mapeamento): %s" % path)
    return data


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes", "on")


@dataclass
class BonderConfig:
    """Estado de configuracao resolvido e pronto para uso."""

    config_dir: Path
    data_root: Path
    runs_dir: Path
    offline: bool
    baits: Dict[str, Any] = field(default_factory=dict)
    samples: Dict[str, Any] = field(default_factory=dict)
    databases: Dict[str, Any] = field(default_factory=dict)
    gpu: Dict[str, Any] = field(default_factory=dict)

    # ---- acessores de conveniencia -----------------------------------------
    @property
    def taxid(self) -> int:
        org = self.databases.get("organism") or self.baits.get("organism_sample") or {}
        return int(org.get("taxid", 10090))

    def bait(self, bait_id: str) -> Dict[str, Any]:
        for b in self.baits.get("baits", []):
            if b.get("id") == bait_id:
                return b
        raise KeyError("Isca desconhecida em baits.yaml: %r" % bait_id)

    def bait_ids(self) -> List[str]:
        return [b.get("id") for b in self.baits.get("baits", [])]

    def contrasts(self) -> List[Dict[str, Any]]:
        return list(self.samples.get("contrasts", []))

    def samples_for(self, condition: str) -> List[Dict[str, Any]]:
        return [s for s in self.samples.get("samples", []) if s.get("condition") == condition]

    def cache_path(self, *parts: str) -> Path:
        p = self.data_root.joinpath(*parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def run_path(self, run_id: str, *parts: str) -> Path:
        p = self.runs_dir.joinpath(run_id, *parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p


def load_config() -> BonderConfig:
    """Carrega e resolve toda a configuracao do bonder."""
    cfg_dir = _config_dir()
    databases = _load_yaml(cfg_dir / "databases.yaml")
    baits = _load_yaml(cfg_dir / "baits.yaml")
    samples = _load_yaml(cfg_dir / "samples.yaml")
    gpu = _load_yaml(cfg_dir / "gpu.yaml") if (cfg_dir / "gpu.yaml").exists() else {}

    # data_root: env > yaml > default; relativo a PROJECT_ROOT quando nao-absoluto.
    dr_raw = os.environ.get("BONDER_DATA_ROOT") or databases.get("data_root") or "data"
    data_root = Path(dr_raw).expanduser()
    if not data_root.is_absolute():
        data_root = (PROJECT_ROOT / data_root).resolve()

    runs_raw = os.environ.get("BONDER_RUNS_DIR") or "runs"
    runs_dir = Path(runs_raw).expanduser()
    if not runs_dir.is_absolute():
        runs_dir = (PROJECT_ROOT / runs_dir).resolve()

    offline = _as_bool(os.environ.get("BONDER_OFFLINE"), _as_bool(databases.get("offline"), False))

    return BonderConfig(
        config_dir=cfg_dir,
        data_root=data_root,
        runs_dir=runs_dir,
        offline=offline,
        baits=baits,
        samples=samples,
        databases=databases,
        gpu=gpu,
    )
