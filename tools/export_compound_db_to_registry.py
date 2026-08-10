# -*- coding: utf-8 -*-
"""
Converte a biblioteca de ligantes em arquivos declarativos de registry/ligands/.

Funde tres fontes, nessa ordem de precedencia:
  1. os metadata.json ja existentes em docking_project/ligand/<nome>/ (dados
     medidos: formula, MW, atomos pesados, torsoes, carga);
  2. compound_db.DEFAULT_COMPOUND (o 2-SCN-NQ, base do SUB01);
  3. compound_db.NQ_LIBRARY (SMILES e nomes de exibicao da serie).

O campo `class` deixa de ser a categoria fixa "naphthoquinone_*" de
populate_ligand_db.py e passa a ser metadado por ligante -- e o que permite a
mesma plataforma abrigar um macrociclo de 19 membros ao lado de uma naftoquinona.

Uso:
    python tools/export_compound_db_to_registry.py --dry-run
    python tools/export_compound_db_to_registry.py
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "registry" / "ligands"
LIGAND_DIRS = [
    ROOT / "docking_project" / "ligand",
]


def registry_id(name):
    """'2-SCN-NQ' -> '2-scn-nq'. Estavel em filesystem case-insensitive."""
    return re.sub(r"[^a-z0-9\-_]+", "-", str(name).strip().lower()).strip("-") or "unnamed"


def load_disk_metadata():
    """metadata.json de cada composto ja preparado em disco, por nome."""
    found = {}
    for base in LIGAND_DIRS:
        if not base.exists():
            continue
        for meta_path in base.glob("*/metadata.json"):
            try:
                with open(meta_path, encoding="utf-8") as fh:
                    meta = json.load(fh)
            except (json.JSONDecodeError, UnicodeDecodeError) as exc:
                print(f"[aviso] {meta_path} ilegivel: {exc}")
                continue
            name = meta.get("name") or meta_path.parent.name
            found[name] = meta
    return found


def build_ligand(name, lib_entry, meta):
    """Monta o ligante declarativo. Campos ausentes ficam None, nunca inventados."""
    smiles = (meta.get("smiles") or (lib_entry or {}).get("smiles") or "").strip()
    rotatable = meta.get("rotatable_bonds")

    payload = {
        "schema_version": "2.0",
        "id": registry_id(name),
        "name": name,
        "display_name": meta.get("display_name") or (lib_entry or {}).get("display_name") or name,
        "smiles": smiles,
        "formula": meta.get("formula"),
        "mw_g_mol": meta.get("mw"),
        "heavy_atoms": meta.get("heavy_atoms"),
        "charge": meta.get("charge", 0),
        "class": "naphthoquinone",
        "chemotype_tags": ["naphthoquinone"],
        "topology": {
            # A serie de NQ nao tem macrociclo; o campo existe para que um
            # ligante que tenha (jasplaquinolida, anel de 19) declare o requisito
            # de flexibilidade em vez de depender do default da ferramenta.
            "is_macrocycle": False,
            "largest_ring_size": None,
            "expected_active_torsions": rotatable,
        },
        "prep": {
            "optimization": meta.get("optimization"),
            "protonation": None,
            "pdbqt": {"tool": "openbabel"},
            "gates": [],
        },
        "notes": meta.get("notes") or (lib_entry or {}).get("notes") or "",
        "provenance": {
            "generated_from": "compound_db + docking_project/ligand/*/metadata.json",
            "generator": "tools/export_compound_db_to_registry.py",
            "date_added": meta.get("date_added"),
            "has_prepared_pdbqt": bool(meta.get("_pdbqt_found")),
        },
    }
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    import compound_db

    disk = load_disk_metadata()
    library = {e["name"]: e for e in compound_db.NQ_LIBRARY}
    default = compound_db.DEFAULT_COMPOUND
    names = sorted(set(library) | set(disk) | {default["name"]})

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written = skipped = 0

    for name in names:
        meta = dict(disk.get(name) or {})
        if name == default["name"]:
            # o default preenche lacunas, sem sobrepor o que foi medido em disco
            for k, v in default.items():
                meta.setdefault(k, v)

        payload = build_ligand(name, library.get(name), meta)
        path = OUT_DIR / f"{payload['id']}.json"

        if path.exists() and not args.force:
            print(f"[skip] {path.name}")
            skipped += 1
            continue
        if args.dry_run:
            print(f"[dry] {path.name}  smiles={payload['smiles'][:40]!r}")
            continue

        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        print(f"[ok] {path.name}")
        written += 1

    missing = [n for n in names if not (disk.get(n, {}).get("smiles") or library.get(n, {}).get("smiles"))]
    if missing:
        print(f"\n[aviso] sem SMILES: {missing}")
    print(f"\n{written} escrito(s), {skipped} pulado(s), {len(names)} ligante(s) no total")


if __name__ == "__main__":
    main()
