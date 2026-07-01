#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cli.py - interface de linha de comando do bonder (headless / CI).

Subcomandos:
  bonder db sync              baixa e verifica os bancos murinos (databases.yaml)
  bonder baits check          valida iscas e monta o FASTA de busca combinado
  bonder qc --input F         resumo de QC de uma saida de busca
  bonder run --input F [...]  roda o pipeline para um contraste (isca vs controle)
  bonder list-contrasts       lista os contrastes definidos em samples.yaml
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from bonder import __version__
from bonder.config import load_config


def _run_id(prefix: str = "run") -> str:
    return "%s_%s" % (prefix, time.strftime("%Y%m%d_%H%M%S"))


def cmd_db_sync(args) -> int:
    from bonder.db.downloader import Downloader
    cfg = load_config()
    dl = Downloader(cfg.data_root, offline=cfg.offline)
    sources = cfg.databases.get("sources", {})
    plan = [
        ("uniprot_fasta", sources.get("uniprot_mouse", {}).get("fasta_url"), "uniprot/mouse.fasta",
         sources.get("uniprot_mouse", {}).get("sha256_fasta")),
        ("uniprot_idmap", sources.get("uniprot_mouse", {}).get("tsv_url"), "uniprot/mouse_idmap.tsv",
         sources.get("uniprot_mouse", {}).get("sha256_tsv")),
        ("string_links", sources.get("string_mouse", {}).get("links_url"),
         "string/10090.links.txt.gz", sources.get("string_mouse", {}).get("sha256_links")),
        ("string_aliases", sources.get("string_mouse", {}).get("aliases_url"),
         "string/10090.aliases.txt.gz", sources.get("string_mouse", {}).get("sha256_aliases")),
        ("contaminants", sources.get("contaminants", {}).get("url"), "contaminants/crap.fasta",
         sources.get("contaminants", {}).get("sha256")),
    ]
    ok = 0
    for name, url, dest, sha in plan:
        if not url:
            print("[skip] %s: url nao definida" % name)
            continue
        try:
            p = dl.fetch(name, url, dest, expected_sha256=sha)
            print("[ok]   %s -> %s" % (name, p))
            ok += 1
        except Exception as e:
            print("[erro] %s: %s" % (name, e))
    print("Sincronizacao: %d/%d fontes prontas em %s" % (ok, len(plan), cfg.data_root))
    return 0 if ok else 1


def cmd_baits_check(args) -> int:
    from bonder.db.baits import validate_baits, build_search_fasta
    cfg = load_config()
    warns = validate_baits(cfg.baits, cfg.config_dir)
    for w in warns:
        print("[aviso] %s" % w)
    fasta = cfg.cache_path("uniprot", "mouse.fasta")
    contam = cfg.cache_path("contaminants", "crap.fasta")
    out = build_search_fasta(
        cfg.baits, cfg.config_dir, cfg.data_root / "search_db.fasta",
        mouse_fasta=fasta if Path(fasta).exists() else None,
        contaminants_fasta=contam if Path(contam).exists() else None,
    )
    print("FASTA de busca: %s" % out)
    return 0


def cmd_qc(args) -> int:
    from bonder.io import load_search_output
    from bonder.stats import qc
    cfg = load_config()
    ds = load_search_output(args.input, samples_cfg=cfg.samples.get("samples", []),
                            baits_cfg=cfg.baits, report_tsv=args.report)
    print(json.dumps(qc.qc_summary(ds), indent=2, ensure_ascii=True))
    return 0


def cmd_list_contrasts(args) -> int:
    cfg = load_config()
    for c in cfg.contrasts():
        print("%-16s bait=%-14s test=%-12s control=%s"
              % (c.get("name"), c.get("bait_id"), c.get("test"), c.get("control")))
    return 0


def cmd_run(args) -> int:
    from bonder.pipeline import run_pipeline
    cfg = load_config()
    contrasts = cfg.contrasts()
    if args.contrast:
        contrasts = [c for c in contrasts if c.get("name") == args.contrast]
    elif args.bait:
        contrasts = [c for c in contrasts if c.get("bait_id") == args.bait]
    if not contrasts:
        print("Nenhum contraste correspondente em samples.yaml.", file=sys.stderr)
        return 2
    for c in contrasts:
        rid = args.run_id or _run_id(c.get("name", "run"))
        print("== rodando contraste %s (run_id=%s) ==" % (c.get("name"), rid))
        run_dir = run_pipeline(args.input, c, rid, cfg=cfg, report_tsv=args.report)
        print("saida: %s" % run_dir)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="bonder", description="Plataforma CoIP-MS / AP-MS (proteoma murino).")
    p.add_argument("--version", action="version", version="bonder %s" % __version__)
    sub = p.add_subparsers(dest="command", required=True)

    db = sub.add_parser("db", help="bancos de referencia")
    dbsub = db.add_subparsers(dest="db_command", required=True)
    dbsub.add_parser("sync", help="baixa e verifica os bancos murinos").set_defaults(func=cmd_db_sync)

    bc = sub.add_parser("baits", help="iscas")
    bcsub = bc.add_subparsers(dest="baits_command", required=True)
    bcsub.add_parser("check", help="valida iscas e monta o FASTA de busca").set_defaults(func=cmd_baits_check)

    q = sub.add_parser("qc", help="resumo de QC")
    q.add_argument("--input", required=True)
    q.add_argument("--report", default=None, help="report.tsv do DIA-NN (opcional)")
    q.set_defaults(func=cmd_qc)

    lc = sub.add_parser("list-contrasts", help="lista contrastes de samples.yaml")
    lc.set_defaults(func=cmd_list_contrasts)

    r = sub.add_parser("run", help="roda o pipeline")
    r.add_argument("--input", required=True, help="saida do motor de busca")
    r.add_argument("--report", default=None, help="report.tsv do DIA-NN (opcional)")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--contrast", default=None, help="nome do contraste (samples.yaml)")
    g.add_argument("--bait", default=None, help="id da isca (samples.yaml)")
    r.add_argument("--run-id", dest="run_id", default=None)
    r.set_defaults(func=cmd_run)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
