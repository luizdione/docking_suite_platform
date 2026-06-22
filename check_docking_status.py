#!/usr/bin/env python3
"""
check_docking_status.py - Verifica disponibilidade de poses docked (result_<TID>.pdbqt)
para cada alvo de MD. Roda os mesmos caminhos de busca que md_module.py usa.

Usage:
    python check_docking_status.py
    python check_docking_status.py --targets TcTR
    python check_docking_status.py --compound 2-SCN-NQ --targets TcTR LiTR hGR
"""

import argparse
import os as _os
from pathlib import Path
import sys


def safe_name(name: str) -> str:
    """Replica compound_db._safe_name."""
    return name.replace("/", "_").replace("\\", "_").replace(" ", "_")


def search_docked_pose(target_id: str, compound_name: str, base_dir: str = "docking_project"):
    """Reproduz exatamente a logica de busca do md_module.run_md_analysis."""
    safe_c = safe_name(compound_name)
    script_dir = Path(__file__).resolve().parent
    base = Path(base_dir)
    base_abs = base.resolve() if base.is_absolute() else (script_dir / base).resolve()

    # NQ_DATA_DIR ou fallback para docking_project relativo ao script
    nq_data = Path(_os.environ.get("NQ_DATA_DIR", "")) if _os.environ.get("NQ_DATA_DIR") else None
    roots = [
        base_abs,
        base.resolve(),
        script_dir / "docking_project",
        script_dir,
    ]
    if nq_data:
        roots += [nq_data, nq_data / "docking_project"]
    seen, roots_unique = set(), []
    for r in roots:
        rk = str(r).lower()
        if rk not in seen:
            seen.add(rk); roots_unique.append(r)

    found = []
    for root in roots_unique:
        if not root.exists():
            continue
        # 1) figures/<compound>/<session>/result_<TID>.pdbqt  (canonical)
        fig_base = root / "figures" / safe_c
        if fig_base.exists():
            for sess in sorted([s for s in fig_base.iterdir() if s.is_dir()], reverse=True):
                cand = sess / f"result_{target_id}.pdbqt"
                if cand.exists():
                    found.append(("docked_session", cand, cand.stat().st_size))
        # 2) figures/result_<TID>.pdbqt (flat fallback)
        flat = root / "figures" / f"result_{target_id}.pdbqt"
        if flat.exists():
            found.append(("docked_flat", flat, flat.stat().st_size))
        # 3) ligand/<compound>/*.pdbqt (preparation - INVALIDO para MD)
        ld = root / "ligand" / safe_c
        if ld.exists():
            for cand in ld.glob("*.pdbqt"):
                found.append(("preparation_ONLY", cand, cand.stat().st_size))
    return found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", nargs="+",
                        default=["TcTR", "LiTR", "TcCruzain", "TcFeSOD",
                                 "LmDHODH", "hGR", "hCatL", "hDHODH", "hMnSOD"])
    parser.add_argument("--compound", default="2-SCN-NQ")
    parser.add_argument("--base-dir", default="docking_project")
    args = parser.parse_args()

    print(f"\n{'='*70}")
    print(f"Diagnostico de poses docked")
    print(f"  Compound: {args.compound}")
    print(f"  Targets:  {', '.join(args.targets)}")
    print(f"{'='*70}\n")

    overall_ok = True
    for tid in args.targets:
        print(f"  [{tid}]")
        results = search_docked_pose(tid, args.compound, args.base_dir)
        if not results:
            print(f"    (nada encontrado em nenhuma raiz de busca)")
            overall_ok = False
            continue
        has_docked = any(s in ("docked_session", "docked_flat") for s, _, _ in results)
        for source, path, size in results:
            tag = "OK" if source.startswith("docked") else "AVISO"
            print(f"    [{tag}] {source:<22} {size:>10} bytes  {path}")
        if not has_docked:
            print(f"    >>> ERRO: Apenas PDBQT de PREPARACAO encontrado.")
            print(f"    >>> Rode primeiro: python docking_v2.py --targets {tid} --exhaustiveness 32")
            overall_ok = False
        print()

    print(f"{'='*70}")
    if overall_ok:
        print(f"OK - Todos os alvos tem pose docked. Pode rodar MD com seguranca.")
        sys.exit(0)
    else:
        # Coletar lista exata dos alvos faltantes para gerar comando pronto
        missing = []
        for tid in args.targets:
            results = search_docked_pose(tid, args.compound, args.base_dir)
            has_docked = any(s in ("docked_session", "docked_flat") for s, _, _ in results)
            if not has_docked:
                missing.append(tid)
        print(f"ATENCAO - {len(missing)} alvo(s) sem pose docked. Rode docking primeiro!")
        if missing:
            tids_str = ",".join(missing)
            print(f"")
            print(f"Comando direto para os alvos faltantes:")
            print(f"  python docking_v2.py --targets {tids_str} --exhaustiveness 32")
            print(f"")
            print(f"Ou rode o batch completo (todos os 9 alvos do estudo):")
            print(f"  .\\run_docking_batch.bat")
        sys.exit(1)


if __name__ == "__main__":
    main()
