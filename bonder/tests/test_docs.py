#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Testes do registro de documentos (aba Ajuda e Tutoriais)."""
from bonder.web import docs


def test_list_docs_has_core_entries():
    ids = {d["id"] for d in docs.list_docs()}
    for expected in ("visao-geral", "como-rodar", "metodologia",
                     "reprodutibilidade", "versionamento", "licenca"):
        assert expected in ids, "faltou doc %s" % expected


def test_get_doc_returns_markdown():
    d = docs.get_doc("metodologia")
    assert d is not None
    assert d["title"] and d["markdown"]
    assert "SAINTexpress" in d["markdown"]  # conteudo real


def test_get_doc_unknown_is_none():
    assert docs.get_doc("nao-existe") is None
    # ids nao vazam caminho: um id arbitrario nao resolve
    assert docs.get_doc("../config/baits") is None


def test_docs_endpoints(tmp_path, monkeypatch):
    from bonder.web import create_app
    c = create_app().test_client()
    r = c.get("/api/docs")
    assert r.status_code == 200
    ids = {d["id"] for d in r.get_json()}
    assert "licenca" in ids and "versionamento" in ids
    rl = c.get("/api/docs/licenca")
    assert rl.status_code == 200 and "MIT" in rl.get_json()["markdown"]
    assert c.get("/api/docs/nao-existe").status_code == 404
