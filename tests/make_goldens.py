# -*- coding: utf-8 -*-
"""
Gera os snapshots golden do SUB01 (Fase 0 da generalizacao multi-alvo).

Estes arquivos congelam o comportamento ATUAL da plataforma, antes de qualquer
refatoracao. Sao o gate de aceite de todas as fases seguintes: o registry
declarativo tem que reproduzi-los byte a byte.

O golden de centros de grade e o mais importante -- o centro da caixa e o numero
que, se mudar, invalida o docking do paper ja submetido (SUB01 / 2-SCN-NQ).

Uso:
    python tests/make_goldens.py            # gera (nao sobrescreve se existir)
    python tests/make_goldens.py --force    # regenera
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"
PROTEIN_DIR = ROOT / "docking_project" / "proteins"


def snapshot_target_db():
    """Serializa TARGET_DB preservando a informacao de tipo das tuplas.

    grid_size e tupla no codigo e vira lista em JSON. Guardamos o tipo
    explicitamente para que o teste de paridade detecte a ida-e-volta.
    """
    import docking_v2

    out = {}
    for tid, info in docking_v2.TARGET_DB.items():
        entry = dict(info)
        entry["_types"] = {k: type(v).__name__ for k, v in info.items()}
        out[tid] = entry
    return {
        "_meta": {
            "source": "docking_v2.TARGET_DB",
            "n_targets": len(out),
            "note": "snapshot pre-refatoracao; nao editar a mao",
        },
        "targets": out,
    }


def snapshot_grid_centers(download_missing=False):
    """Roda _calculate_grid_center() sobre os PDBs versionados no repo."""
    import docking_v2
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    centers, missing = {}, []

    for tid, info in docking_v2.TARGET_DB.items():
        pdb_path = PROTEIN_DIR / f"{info['pdb']}.pdb"

        if not pdb_path.exists() and download_missing:
            got, msg = docking_v2.download_pdb(info["pdb"], PROTEIN_DIR)
            print(f"  download {info['pdb']}: {msg}")
            if got:
                pdb_path = got

        if not pdb_path.exists():
            missing.append({"target": tid, "pdb": info["pdb"]})
            continue

        structure = parser.get_structure("prot", str(pdb_path))
        center = docking_v2._calculate_grid_center(structure, info)
        centers[tid] = {
            "pdb": info["pdb"],
            "grid_strategy": info["grid_strategy"],
            "grid_center": list(center),
            "grid_size": list(info["grid_size"]),
            "pdb_sha256": _sha256(pdb_path),
        }

    return {
        "_meta": {
            "source": "docking_v2._calculate_grid_center",
            "protein_dir": str(PROTEIN_DIR.relative_to(ROOT)).replace("\\", "/"),
            "n_covered": len(centers),
            "missing_pdb_locally": missing,
            "note": (
                "GATE BLOQUEANTE. Se um destes centros mudar, o docking do SUB01 "
                "deixa de ser reproduzivel e o paper submetido fica orfao."
            ),
        },
        "centers": centers,
    }


def _sha256(path):
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write(path, payload, force):
    if path.exists() and not force:
        print(f"[skip] {path.name} ja existe (use --force para regenerar)")
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, sort_keys=True)
        fh.write("\n")
    print(f"[ok] {path.name}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument(
        "--download-missing",
        action="store_true",
        help="baixa do RCSB os PDBs de alvos ausentes localmente",
    )
    args = ap.parse_args()

    print("Gerando goldens do SUB01...")
    write(GOLDEN_DIR / "sub01_targets_golden.json", snapshot_target_db(), args.force)
    write(
        GOLDEN_DIR / "sub01_grid_centers_golden.json",
        snapshot_grid_centers(args.download_missing),
        args.force,
    )


if __name__ == "__main__":
    main()
