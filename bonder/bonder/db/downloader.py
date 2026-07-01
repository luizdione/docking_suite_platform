#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""downloader.py - download reproduzivel e verificado dos bancos de referencia.

Principios (mesma filosofia _deps da suite: nada silencioso, tudo reproduzivel):
  - versoes FIXADAS em config/databases.yaml (nunca 'latest');
  - escrita atomica (tmp -> rename) e retry com backoff exponencial;
  - registry.json no cache com url, sha256, tamanho e timestamp;
  - checksum esperado (quando definido) e verificado; divergencia -> erro;
  - modo offline (BONDER_OFFLINE=1): usa so o cache; erro claro se faltar.

Rode 'bonder db sync' uma vez para popular o cache antes das analises.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import requests
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["requests"])
    import requests


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class Downloader:
    """Cache content-addressed com registry e verificacao de integridade."""

    def __init__(self, data_root: Path, offline: bool = False):
        self.root = Path(data_root)
        self.offline = offline
        self.root.mkdir(parents=True, exist_ok=True)
        self.registry_path = self.root / "registry.json"
        self.registry: Dict[str, Any] = {}
        if self.registry_path.exists():
            try:
                self.registry = json.loads(self.registry_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                self.registry = {}

    def _save_registry(self) -> None:
        self.registry_path.write_text(
            json.dumps(self.registry, indent=2, ensure_ascii=True), encoding="utf-8"
        )

    def fetch(
        self,
        name: str,
        url: str,
        dest_rel: str,
        expected_sha256: Optional[str] = None,
        retries: int = 4,
    ) -> Path:
        """Baixa `url` para `dest_rel` (relativo ao cache) de forma verificada.

        Se o arquivo ja existe e o sha256 confere (ou nao ha esperado), reusa.
        Em modo offline, exige que ja esteja em cache.
        """
        dest = self.root / dest_rel
        dest.parent.mkdir(parents=True, exist_ok=True)

        if dest.exists():
            got = sha256_of(dest)
            if expected_sha256 and got != expected_sha256:
                raise ValueError(
                    "Checksum divergente para %s em cache (esperado %s, obtido %s). "
                    "Remova o arquivo e re-sincronize." % (name, expected_sha256, got)
                )
            self._record(name, url, dest, got)
            return dest

        if self.offline:
            raise FileNotFoundError(
                "Modo offline: %s ausente do cache (%s). Rode 'bonder db sync' "
                "com rede habilitada antes de usar offline." % (name, dest)
            )

        tmp = dest.with_suffix(dest.suffix + ".partial")
        last_err: Optional[Exception] = None
        for attempt in range(retries):
            try:
                with requests.get(url, stream=True, timeout=120) as r:
                    r.raise_for_status()
                    with open(tmp, "wb") as fh:
                        for block in r.iter_content(chunk_size=1 << 20):
                            if block:
                                fh.write(block)
                got = sha256_of(tmp)
                if expected_sha256 and got != expected_sha256:
                    tmp.unlink(missing_ok=True)
                    raise ValueError(
                        "Checksum divergente ao baixar %s (esperado %s, obtido %s)."
                        % (name, expected_sha256, got)
                    )
                tmp.replace(dest)  # rename atomico
                self._record(name, url, dest, got)
                return dest
            except (requests.RequestException, ValueError) as e:
                last_err = e
                tmp.unlink(missing_ok=True)
                if isinstance(e, ValueError):
                    raise  # checksum: nao adianta repetir
                time.sleep(2 ** attempt)  # backoff 1,2,4,8s
        raise RuntimeError("Falha ao baixar %s apos %d tentativas: %s" % (name, retries, last_err))

    def _record(self, name: str, url: str, dest: Path, sha: str) -> None:
        self.registry[name] = {
            "url": url,
            "path": str(dest.relative_to(self.root)),
            "sha256": sha,
            "size": dest.stat().st_size,
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        self._save_registry()
