# -*- coding: utf-8 -*-
"""
Auditor de hard-code de dominio (Fase 0 da generalizacao multi-alvo).

Conta, arquivo a arquivo, as suposicoes soldadas no codigo que impedem a
plataforma de ser guarda-chuva de alvos quaisquer: nome do composto do SUB01,
alvos do painel, classe quimica e organismo.

E somente leitura. Serve de checklist vivo: a contagem tem que cair
monotonicamente ao longo das fases. Grave uma linha de base com --save e compare
depois com --check.

Uso:
    python tools/audit_hardcode.py                 # relatorio
    python tools/audit_hardcode.py --save          # grava a linha de base
    python tools/audit_hardcode.py --check         # falha se a divida subiu
    python tools/audit_hardcode.py --detail        # mostra arquivo:linha
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASELINE = Path(__file__).resolve().parent / "hardcode_baseline.json"

# Arquivos onde o hard-code de dominio importa. Os demais modulos
# (admet, toxicity, drug_design...) sao agnosticos a alvo.
WATCHED = [
    "dashboard_v2.py",
    # A UI saiu de dashboard_v2.py para ca na extracao da Fase 5.1. Precisa
    # continuar sendo contada: mudar de arquivo nao paga divida nenhuma, e uma
    # metrica que caisse por causa da mudanca estaria mentindo.
    "static/app.html",
    "docking_v2.py",
    "md_simulation.py",
    "target_selector.py",
    "compound_db.py",
    "populate_ligand_db.py",
    "export_web.py",
    "analysis.py",
]

# categoria -> regex. Nomeadas pela razao de existirem, nao pelo padrao.
PATTERNS = {
    "composto_sub01": re.compile(r"2-SCN-NQ|2_SCN_NQ|2-Thiocyanate", re.I),
    "alvos_do_painel": re.compile(
        r"\b(TcTR|LiTR|hGR|TcCruzain|TcFeSOD|LmPTR1|LmDHODH|hCatL|hDHODH|hMnSOD)\b"
    ),
    "classe_quimica": re.compile(r"naphthoquinone|naftoquinona|\bNQ_LIBRARY\b", re.I),
    "organismo": re.compile(r"T\.\s*cruzi|Trypanosoma|Leishmania|tripanossom", re.I),
    "cofator_flavina": re.compile(r"\bFAD\b|\bFMN\b"),
}


def scan(detail=False):
    counts, hits = {}, {}
    for rel in WATCHED:
        path = ROOT / rel
        if not path.exists():
            counts[rel] = {"_status": "ausente"}
            continue
        per_file = {}
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except UnicodeDecodeError:
            lines = path.read_text(encoding="cp1252", errors="replace").splitlines()
        for cat, rx in PATTERNS.items():
            found = [(i + 1, ln.strip()[:110]) for i, ln in enumerate(lines) if rx.search(ln)]
            per_file[cat] = len(found)
            if detail and found:
                hits.setdefault(rel, {})[cat] = found
        counts[rel] = per_file
    return counts, hits


def total(counts):
    return sum(
        v for f in counts.values() for k, v in f.items() if not k.startswith("_") and isinstance(v, int)
    )


def report(counts, hits, detail):
    cats = list(PATTERNS)
    width = max(len(f) for f in counts) + 2
    header = "arquivo".ljust(width) + "".join(c[:14].rjust(16) for c in cats) + "total".rjust(8)
    print(header)
    print("-" * len(header))
    for rel, per in sorted(counts.items(), key=lambda kv: -total({kv[0]: kv[1]})):
        if per.get("_status") == "ausente":
            print(rel.ljust(width) + "  (ausente)")
            continue
        row = rel.ljust(width) + "".join(str(per.get(c, 0)).rjust(16) for c in cats)
        print(row + str(sum(per.values())).rjust(8))
    print("-" * len(header))
    print("TOTAL".ljust(width) + "".join(
        str(sum(p.get(c, 0) for p in counts.values() if "_status" not in p)).rjust(16) for c in cats
    ) + str(total(counts)).rjust(8))

    if detail:
        for rel, per_cat in hits.items():
            print(f"\n=== {rel} ===")
            for cat, found in per_cat.items():
                print(f"  [{cat}]")
                for lineno, text in found[:40]:
                    print(f"    {rel}:{lineno}: {text}")
                if len(found) > 40:
                    print(f"    ... +{len(found) - 40} ocorrencias")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", action="store_true", help="grava a linha de base")
    ap.add_argument("--check", action="store_true", help="falha se a divida subiu")
    ap.add_argument("--detail", action="store_true", help="mostra arquivo:linha")
    args = ap.parse_args()

    counts, hits = scan(detail=args.detail)
    report(counts, hits, args.detail)
    now = total(counts)

    if args.save:
        BASELINE.write_text(
            json.dumps({"total": now, "per_file": counts}, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"\n[ok] linha de base gravada: {now} ocorrencias")
        return 0

    if args.check:
        if not BASELINE.exists():
            print("\n[erro] sem linha de base -- rode --save primeiro")
            return 2
        base = json.loads(BASELINE.read_text(encoding="utf-8"))["total"]
        delta = now - base
        print(f"\nlinha de base: {base} | agora: {now} | delta: {delta:+d}")
        if delta > 0:
            print("[FALHA] a divida de hard-code SUBIU")
            return 1
        print("[ok] divida estavel ou em queda")
    return 0


if __name__ == "__main__":
    sys.exit(main())
