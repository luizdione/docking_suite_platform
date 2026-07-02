#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Testes do modulo de manutencao/saude e dos endpoints de auditoria."""
from bonder.config import PROJECT_ROOT


def _env(monkeypatch, tmp_path):
    monkeypatch.setenv("BONDER_CONFIG_DIR", str(PROJECT_ROOT / "config"))
    monkeypatch.setenv("BONDER_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("BONDER_DATA_ROOT", str(tmp_path / "data"))
    monkeypatch.delenv("BONDER_REDIS_URL", raising=False)


def test_full_report_structure(tmp_path, monkeypatch):
    _env(monkeypatch, tmp_path)
    from bonder.web import health
    r = health.full_report()
    for key in ("overall", "system", "python_deps", "tools", "r_packages",
                "gpu", "config", "databases", "paths", "runs"):
        assert key in r, "faltou secao %s" % key
    # nucleo (numpy/pandas/pyarrow/networkx/yaml/requests/flask) presente no ambiente de teste
    assert r["overall"] == "ok", r["core_missing"]
    assert "NSs_vs_CTRL" in r["config"]["contrasts"]
    # sem GPU no CI
    assert r["gpu"]["device"]["device"] == "cpu"
    # ferramentas externas ausentes no CI, mas listadas
    assert any(t["name"] == "SAINTexpress-int" for t in r["tools"])


def test_self_test_passes(tmp_path, monkeypatch):
    _env(monkeypatch, tmp_path)
    from bonder.web import health
    s = health.self_test()
    assert s["ok"] is True, s.get("error")
    assert s["n_candidates"] >= 1
    steps = {x["step"]: x["status"] for x in s["steps"]}
    assert steps["ingest"] == "ok"
    assert steps["saint"] in ("pulado", "erro")  # binario ausente no CI
    # bait e contaminante nao entram como candidatos
    arts = {c["artifact"]: c["ok"] for c in s["checks"]}
    assert arts["bait excluida dos preys"] and arts["contaminante excluido"]


def test_health_endpoints(tmp_path, monkeypatch):
    _env(monkeypatch, tmp_path)
    from bonder.web import create_app
    c = create_app().test_client()
    rf = c.get("/api/health/full")
    assert rf.status_code == 200 and rf.get_json()["overall"] == "ok"
    rs = c.post("/api/health/selftest")
    assert rs.status_code == 200 and rs.get_json()["ok"] is True
