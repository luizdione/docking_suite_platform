#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Testa o caminho de job assincrono (backend thread, sem Redis)."""
import time
from pathlib import Path

from bonder.config import PROJECT_ROOT
from tests.test_pipeline_e2e import _write_maxquant_full


def test_thread_job_runs(tmp_path, monkeypatch):
    monkeypatch.setenv("BONDER_CONFIG_DIR", str(PROJECT_ROOT / "config"))
    monkeypatch.setenv("BONDER_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("BONDER_DATA_ROOT", str(tmp_path / "data"))
    monkeypatch.delenv("BONDER_REDIS_URL", raising=False)

    from bonder.config import load_config
    from bonder.web import jobs
    cfg = load_config()
    mq = tmp_path / "proteinGroups.txt"
    _write_maxquant_full(mq)
    contrast = next(c for c in cfg.contrasts() if c["name"] == "NSs_vs_CTRL")

    backend = jobs.enqueue_run("web_run", str(mq), contrast)
    assert backend == "thread"

    for _ in range(100):  # ate ~10s
        st = jobs.read_status("web_run")
        if st and st["status"] in ("done", "error"):
            break
        time.sleep(0.1)
    st = jobs.read_status("web_run")
    assert st["status"] == "done", st
    assert (cfg.runs_dir / "web_run" / "interactors.tsv").exists()
