#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""app.py - dashboard Flask do bonder (localhost).

Segue a convencao Flask da suite docking (dashboard_v2). Sobe uma UI de navegador
e uma API REST para: carregar uma saida de busca, escolher o contraste, disparar
a corrida (job assincrono) e visualizar QC, tabela de interatores, volcano e a
rede de interacao. Uso local mono-usuario (bind em 127.0.0.1 por padrao).
"""
from __future__ import annotations

import time
from pathlib import Path

try:
    from flask import Flask, jsonify, render_template, request, send_from_directory, abort
except ImportError:  # pragma: no cover
    from bonder import _deps
    _deps.ensure(["flask"])
    from flask import Flask, jsonify, render_template, request, send_from_directory, abort

from bonder import __version__
from bonder.config import load_config
from bonder.web import jobs


def create_app() -> "Flask":
    app = Flask(__name__, template_folder="templates", static_folder="static")

    @app.get("/")
    def index():
        cfg = load_config()
        return render_template("dashboard.html", version=__version__,
                               contrasts=cfg.contrasts(), baits=cfg.baits.get("baits", []))

    @app.get("/api/health")
    def health():
        return jsonify({"status": "ok", "version": __version__, "time": time.strftime("%Y-%m-%dT%H:%M:%S")})

    @app.get("/api/contrasts")
    def contrasts():
        return jsonify(load_config().contrasts())

    @app.get("/api/gpu")
    def gpu():
        from bonder.gpu.detect import probe
        from bonder.gpu import gpu_gate
        return jsonify({"device": probe().to_dict(), "gpu_free_mb": gpu_gate.gpu_free_mb(),
                        "llm_present": gpu_gate.llm_present()})

    @app.post("/api/run")
    def run():
        data = request.get_json(force=True, silent=True) or {}
        input_path = data.get("input")
        contrast_name = data.get("contrast")
        if not input_path or not Path(input_path).exists():
            return jsonify({"error": "input inexistente: %r" % input_path}), 400
        cfg = load_config()
        contrast = next((c for c in cfg.contrasts() if c.get("name") == contrast_name), None)
        if contrast is None:
            return jsonify({"error": "contraste desconhecido: %r" % contrast_name}), 400
        run_id = "run_%s_%s" % (contrast_name, time.strftime("%Y%m%d_%H%M%S"))
        backend = jobs.enqueue_run(run_id, input_path, contrast, report_tsv=data.get("report"))
        return jsonify({"run_id": run_id, "backend": backend, "status": "queued"})

    @app.get("/api/runs")
    def list_runs():
        cfg = load_config()
        runs = []
        if cfg.runs_dir.exists():
            for d in sorted(cfg.runs_dir.iterdir(), reverse=True):
                if d.is_dir():
                    st = jobs.read_status(d.name) or {}
                    runs.append({"run_id": d.name, "status": st.get("status", "unknown"),
                                 "contrast": st.get("contrast")})
        return jsonify(runs)

    @app.get("/api/runs/<run_id>")
    def run_detail(run_id):
        cfg = load_config()
        rd = cfg.runs_dir / run_id
        if not rd.exists():
            abort(404)
        status = jobs.read_status(run_id) or {"status": "unknown"}
        results = {}
        rj = rd / "results.json"
        if rj.exists():
            import json
            results = json.loads(rj.read_text(encoding="utf-8"))
        return jsonify({"status": status, "results": results})

    @app.get("/api/runs/<run_id>/network")
    def run_network(run_id):
        cfg = load_config()
        f = cfg.runs_dir / run_id / "network.json"
        if not f.exists():
            abort(404)
        return send_from_directory(f.parent, f.name, mimetype="application/json")

    @app.get("/api/runs/<run_id>/volcano.png")
    def run_volcano(run_id):
        cfg = load_config()
        f = cfg.runs_dir / run_id / "volcano.png"
        if not f.exists():
            abort(404)
        return send_from_directory(f.parent, f.name, mimetype="image/png")

    @app.get("/api/runs/<run_id>/interactors")
    def run_interactors(run_id):
        cfg = load_config()
        f = cfg.runs_dir / run_id / "interactors.tsv"
        if not f.exists():
            abort(404)
        return send_from_directory(f.parent, f.name, mimetype="text/tab-separated-values")

    return app


def main():
    import os
    app = create_app()
    host = os.environ.get("BONDER_HOST", "127.0.0.1")
    port = int(os.environ.get("BONDER_PORT", "8080"))
    app.run(host=host, port=port, debug=bool(os.environ.get("BONDER_DEBUG")))


if __name__ == "__main__":
    main()
