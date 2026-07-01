#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""jobs.py - execucao assincrona das corridas do pipeline.

Dois backends, escolhidos automaticamente:
  - RQ + Redis (BONDER_REDIS_URL definido e rq/redis disponiveis): fila real,
    consumida pelo conteiner 'worker' (recomendado no docker compose);
  - Thread local (fallback): roda em uma thread no proprio processo Flask, util
    para uso mono-usuario em localhost sem Redis.

O status de cada corrida e persistido em <run_dir>/status.json para o dashboard.
"""
from __future__ import annotations

import json
import os
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Dict, Optional

from bonder.config import load_config


def _status_path(run_id: str) -> Path:
    cfg = load_config()
    return cfg.run_path(run_id, "status.json")


def write_status(run_id: str, status: str, **extra) -> None:
    p = _status_path(run_id)
    payload = {"run_id": run_id, "status": status, "updated": time.strftime("%Y-%m-%dT%H:%M:%S"), **extra}
    p.write_text(json.dumps(payload, indent=2, ensure_ascii=True, default=str), encoding="utf-8")


def read_status(run_id: str) -> Optional[Dict[str, Any]]:
    p = _status_path(run_id)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def execute_run(run_id: str, input_path: str, contrast: Dict[str, Any], report_tsv: Optional[str] = None) -> str:
    """Executa o pipeline e atualiza o status. Usado por RQ ou pela thread local."""
    from bonder.pipeline import run_pipeline
    write_status(run_id, "running", contrast=contrast.get("name"))
    try:
        cfg = load_config()
        run_dir = run_pipeline(input_path, contrast, run_id, cfg=cfg, report_tsv=report_tsv)
        write_status(run_id, "done", contrast=contrast.get("name"), run_dir=str(run_dir))
        return str(run_dir)
    except Exception as e:  # nunca deixa a corrida em estado ambiguo
        write_status(run_id, "error", contrast=contrast.get("name"),
                     error=str(e), traceback=traceback.format_exc())
        raise


def _rq_queue():
    url = os.environ.get("BONDER_REDIS_URL")
    if not url:
        return None
    try:
        from redis import Redis
        from rq import Queue
        return Queue("bonder", connection=Redis.from_url(url))
    except Exception:
        return None


def enqueue_run(run_id: str, input_path: str, contrast: Dict[str, Any], report_tsv: Optional[str] = None) -> str:
    """Enfileira (RQ) ou dispara (thread) uma corrida. Retorna o backend usado."""
    write_status(run_id, "queued", contrast=contrast.get("name"))
    q = _rq_queue()
    if q is not None:
        q.enqueue(execute_run, run_id, input_path, contrast, report_tsv, job_timeout=3600)
        return "rq"
    t = threading.Thread(target=execute_run, args=(run_id, input_path, contrast, report_tsv), daemon=True)
    t.start()
    return "thread"
