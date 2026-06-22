#!/usr/bin/env python3
"""
=============================================================================
COMPOUND DATABASE MODULE
Manages the ligand library in docking_project/ligand/
Each compound lives in its own subdirectory:
  docking_project/ligand/{compound_name}/
    metadata.json   — name, smiles, formula, mw, heavy_atoms, ...
    {name}.pdbqt    — prepared ligand for Vina
    {name}.mol2     — 3D structure
    {name}.sdf      — 3D structure (optional)
    {name}.svg      — 2D structure image (optional)
=============================================================================
"""
import json
import math
import sys
from pathlib import Path
from datetime import datetime
import logging

log = logging.getLogger(__name__)

import os as _os
# BASE_DIR responde ao mesmo env var que o dashboard usa (NQ_DATA_DIR).
# Fallback: docking_project relativo ao script (compatibilidade legada).
_env_data = _os.environ.get("NQ_DATA_DIR", "")
BASE_DIR   = Path(_env_data) if _env_data else Path(__file__).resolve().parent / "docking_project"
LIGAND_DIR = BASE_DIR / "ligand"

# ─── Biblioteca built-in de naftoquinonas antiparasitárias ──────
NQ_LIBRARY = [
    {"name": "2-SCN-NQ",    "smiles": "O=C1C=C(SC#N)C(=O)c2ccccc21",           "display_name": "2-Tiocianato-1,4-naftoquinona",   "notes": "Composto base do projeto; substrato subversivo/inibidor TR."},
    {"name": "2-Cl-NQ",     "smiles": "O=C1C=C(Cl)C(=O)c2ccccc21",             "display_name": "2-Cloro-1,4-naftoquinona",         "notes": "De Castro et al. IC50 ~2 µM vs T. cruzi."},
    {"name": "2-Br-NQ",     "smiles": "O=C1C=C(Br)C(=O)c2ccccc21",             "display_name": "2-Bromo-1,4-naftoquinona",         "notes": "Eletrofílico; IC50 ~0,8–2 µM vs T. cruzi (De Castro et al.)."},
    {"name": "2-NH2-NQ",    "smiles": "O=C1C=C(N)C(=O)c2ccccc21",              "display_name": "2-Amino-1,4-naftoquinona",         "notes": "Doa H-bond para Glu18 de TcTR; série bem estudada."},
    {"name": "2-OH-NQ",     "smiles": "O=C1C=C(O)C(=O)c2ccccc21",              "display_name": "Lawsona (2-OH-1,4-NQ)",            "notes": "Produto natural; atividade antiparasitária documentada."},
    {"name": "2-Me-NQ",     "smiles": "O=C1C=C(C)C(=O)c2ccccc21",              "display_name": "Menadiona (Vitamina K3)",           "notes": "2-Metil-NQ; substrato subversivo clássico de TR."},
    {"name": "2-NO2-NQ",    "smiles": "O=C1C=C([N+](=O)[O-])C(=O)c2ccccc21",   "display_name": "2-Nitro-1,4-naftoquinona",         "notes": "Grupo saída para conversão a amino."},
    {"name": "2-OMe-NQ",    "smiles": "O=C1C=C(OC)C(=O)c2ccccc21",             "display_name": "2-Metóxi-1,4-naftoquinona",        "notes": "Lopes et al. 1978; ativa vs T. cruzi epimastigotas."},
    {"name": "2-triazol-NQ","smiles": "O=C1C=C(n2cnnn2)C(=O)c2ccccc21",        "display_name": "2-(1,2,4-Triazol-1-il)-NQ",       "notes": "Bahia et al. 2016 MedChemComm; IC50 0,8–3,1 µM T. cruzi."},
    {"name": "2-imidazol-NQ","smiles": "O=C1C=C(n2ccnc2)C(=O)c2ccccc21",       "display_name": "2-(Imidazol-1-il)-NQ",             "notes": "SNAr de 2-halo-NQ + imidazol; naftoimidazol-like."},
    {"name": "beta-lapachone","smiles": "O=C1CC(=O)c2cc3ccc(OC3(C)C)cc2C1",     "display_name": "Beta-lapachona",                   "notes": "Produto natural do lapachol; substrato subversivo de NQO1 e TR."},
    {"name": "plumbagin",    "smiles": "O=C1c2c(O)cccc2C(=O)C=C1C",            "display_name": "Plumbagina",                       "notes": "5-Hidroxi-2-metil-NQ; atividade antiparasitária e antitumoral."},
    {"name": "NQ-SPM",       "smiles": "O=C1C=C(SCCCNCCCCNCCCN)C(=O)c2ccccc21","display_name": "NQ-Espermidina (design dual)",     "notes": "Dual design: NQ + cauda espermidina; substrato subversivo + inibidor TR."},
    {"name": "NQ-ACR",       "smiles": "O=C1C=C(NCCCNc3c4ccccc4nc4ccccc34)C(=O)c2ccccc21","display_name": "NQ-Acridina (mepacrina-mimétic)","notes": "Dual design: NQ + acridina; substrato subversivo + inibidor TR sítio Z."},
    {"name": "BIS-NQ",       "smiles": "O=C1C=C(SCCNCCCCNCCSC3=CC(=O)c4ccccc4C3=O)C(=O)c2ccccc21","display_name": "Bis-naftoquinona (bivalente)","notes": "Dual design: dois núcleos NQ; cruza os dois sub-sítios da TR."},
]

# ─── Default 2-SCN-NQ metadata (from ligand_metadata.json) ───────
DEFAULT_COMPOUND = {
    "name":           "2-SCN-NQ",
    "display_name":   "2-Thiocyanato-1,4-naphthoquinone",
    "smiles":         "O=C1C=C(SC#N)C(=O)c2ccccc21",
    "formula":        "C11H5NO2S",
    "mw":             215.23,
    "heavy_atoms":    15,
    "rotatable_bonds": 1,
    "charge":         0,
    "optimization":   "MMFF94 (1000 steps conjugate gradients)",
    "date_added":     "2024-01-01",
    "notes":          "Base compound of the Santos et al. antiparasitic study.",
}


# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def get_ligand_dir(compound_name: str) -> Path:
    """Return the directory for a given compound."""
    return LIGAND_DIR / _safe_name(compound_name)


def _safe_name(name) -> str:
    """Nome de diretorio seguro para o filesystem.

    Rede de defesa contra entradas invalidas vindas do frontend: objetos JS
    serializados chegavam como '[object Object]' e viravam a pasta espuria
    '_object_Object_'. Rejeita nao-strings, vazios e sentinelas de serializacao,
    devolvendo 'compound' como fallback neutro. Nomes validos (ex.: '2-SCN-NQ')
    passam inalterados.
    """
    import re
    if not isinstance(name, str):
        return "compound"
    s = name.strip()
    if not s or s.lower() in ("[object object]", "undefined", "null", "nan", "none"):
        return "compound"
    safe = re.sub(r'[^\w\-.]', '_', s)
    # rejeita resultados degenerados (so '_'/'.' ou o classico _object_Object_)
    if not safe.strip("_.") or safe.lower().startswith("_object_object"):
        return "compound"
    return safe


def list_compounds() -> list[dict]:
    """
    Scan docking_project/ligand/ and return all compounds with their metadata.
    Returns list of metadata dicts, sorted by date_added then name.
    """
    compounds = []
    if not LIGAND_DIR.exists():
        return compounds
    for subdir in sorted(LIGAND_DIR.iterdir()):
        if not subdir.is_dir():
            continue
        meta_path = subdir / "metadata.json"
        if meta_path.exists():
            try:
                with open(meta_path, encoding='utf-8') as f:
                    meta = json.load(f)
                meta["_dir"] = str(subdir)
                meta["_has_pdbqt"] = any(subdir.glob("*.pdbqt"))
                meta["_has_mol2"]  = any(subdir.glob("*.mol2"))
                compounds.append(meta)
            except Exception as e:
                log.warning("list_compounds: metadata ilegivel em %s: %s", meta_path, e)
        else:
            # Directory exists but no metadata — infer minimal info
            compounds.append({
                "name":    subdir.name,
                "smiles":  "",
                "formula": "",
                "mw":      0,
                "_dir":    str(subdir),
                "_has_pdbqt": any(subdir.glob("*.pdbqt")),
                "_has_mol2":  any(subdir.glob("*.mol2")),
            })
    return compounds


def load_compound(name: str) -> dict | None:
    """Load metadata for a specific compound. Returns None if not found."""
    cdir = get_ligand_dir(name)
    meta_path = cdir / "metadata.json"
    if not meta_path.exists():
        # Try scanning all compounds by name field
        for c in list_compounds():
            if c.get("name", "").lower() == name.lower():
                return c
        return None
    with open(meta_path, encoding='utf-8') as f:
        meta = json.load(f)
    meta["_dir"] = str(cdir)
    meta["_has_pdbqt"] = any(cdir.glob("*.pdbqt"))
    meta["_has_mol2"]  = any(cdir.glob("*.mol2"))
    return meta


def get_pdbqt_path(compound_name: str) -> Path | None:
    """Return path to the .pdbqt file for a compound, or None if missing."""
    cdir = get_ligand_dir(compound_name)
    pdbqt_files = list(cdir.glob("*.pdbqt"))
    if pdbqt_files:
        return pdbqt_files[0]
    return None


def add_compound(name: str, smiles: str, formula: str = "",
                 display_name: str = "", mw: float = 0.0,
                 notes: str = "", prepare_3d: bool = True) -> dict:
    """
    Add a new compound to the database.
    Creates the directory, saves metadata, and optionally generates 3D structure.
    Returns the metadata dict with status info.
    """
    safe = _safe_name(name)
    cdir = LIGAND_DIR / safe
    cdir.mkdir(parents=True, exist_ok=True)

    # Calculate MW from SMILES if not provided
    if mw == 0.0 or formula == "":
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if mw == 0.0:
                    mw = round(Descriptors.MolWt(mol), 3)
                if formula == "":
                    formula = rdMolDescriptors.CalcMolFormula(mol)
        except ImportError:
            log.debug("add_compound: RDKit ausente; MW/formula nao calculados de SMILES")

    meta = {
        "name":            name,
        "display_name":    display_name or name,
        "smiles":          smiles,
        "formula":         formula,
        "mw":              mw,
        "heavy_atoms":     _count_heavy_atoms(smiles),
        "rotatable_bonds": _count_rotatable_bonds(smiles),
        "charge":          0,
        "optimization":    "MMFF94 (1000 steps conjugate gradients)",
        "date_added":      datetime.now().strftime("%Y-%m-%d"),
        "notes":           notes,
    }

    # Save metadata
    with open(cdir / "metadata.json", "w", encoding='utf-8') as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    # Prepare 3D structure if requested
    status_msgs = []
    if prepare_3d:
        pdbqt_path = cdir / f"{safe}.pdbqt"
        mol2_path  = cdir / f"{safe}.mol2"
        if not pdbqt_path.exists():
            ok, msg = generate_pdbqt(smiles, pdbqt_path, mol2_path)
            status_msgs.append(msg)
        else:
            status_msgs.append(f"[OK] PDBQT already exists: {pdbqt_path.name}")

    # Generate 2D SVG structure image
    svg_path = cdir / f"{safe}.svg"
    if not svg_path.exists():
        _generate_svg(smiles, svg_path)

    meta["_dir"]        = str(cdir)
    meta["_has_pdbqt"]  = any(cdir.glob("*.pdbqt"))
    meta["_has_mol2"]   = any(cdir.glob("*.mol2"))
    meta["_status_msgs"] = status_msgs
    return meta


def generate_pdbqt(smiles: str, pdbqt_path: Path, mol2_path: Path = None) -> tuple[bool, str]:
    """
    Generate 3D structure (RDKit ETKDG+MMFF) and convert to PDBQT (obabel CLI).
    Falls back to pybel if obabel not found.  Returns (success, message).
    """
    # ── Método preferido: RDKit 3D + obabel CLI (robusto no env datascience) ──
    try:
        import subprocess, tempfile
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("SMILES inválido")
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception as e:
            log.debug("generate_pdbqt: MMFF optimize falhou (%s); seguindo sem minimizar", e)
        pdb_tmp = str(pdbqt_path).replace(".pdbqt", "_3d.pdb")
        Chem.MolToPDBFile(mol, pdb_tmp)
        # obabel: PDB -> pdbqt
        obabel_exe = _find_obabel()
        r = subprocess.run([obabel_exe, pdb_tmp, "-opdbqt", "-O", str(pdbqt_path)],
                           capture_output=True, text=True, timeout=120)
        if Path(pdbqt_path).exists() and Path(pdbqt_path).stat().st_size > 0:
            if mol2_path:
                subprocess.run([obabel_exe, pdb_tmp, "-omol2", "-O", str(mol2_path)],
                               capture_output=True, text=True, timeout=60)
            return True, f"[OK] {pdbqt_path.name} (RDKit+obabel)"
        return False, f"[ERRO] obabel não gerou pdbqt: {r.stderr[:200]}"
    except Exception:
        pass  # fallback para pybel abaixo (RDKit/obabel indisponivel)

    # ── Fallback: pybel (pode falhar no conda se BABEL_LIBDIR não setado) ──
    try:
        from openbabel import pybel
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D(forcefield="mmff94", steps=500)

        # Energy minimisation
        ff = pybel._forcefields.get("mmff94")
        if ff and ff.Setup(mol.OBMol):
            ff.ConjugateGradients(1000, 1.0e-6)
            ff.GetCoordinates(mol.OBMol)

        mol.write("pdbqt", str(pdbqt_path), overwrite=True)
        if mol2_path:
            mol.write("mol2", str(mol2_path), overwrite=True)

        n_atoms = mol.OBMol.NumHvyAtoms()
        size    = pdbqt_path.stat().st_size if pdbqt_path.exists() else 0
        return True, f"[OK] 3D generated: {n_atoms} heavy atoms, {size:,} bytes"

    except ImportError:
        return False, "[ERRO] OpenBabel not installed (pip install openbabel-wheel)"
    except Exception as e:
        return False, f"[ERRO] 3D generation failed: {e}"


def _find_obabel() -> str:
    """Retorna o caminho do executável obabel (datascience env ou PATH)."""
    import subprocess
    candidates = [
        r"C:\Users\luizd\miniforge3\envs\datascience\Library\bin\obabel.exe",
        r"C:\Users\luizd\miniforge3\Library\bin\obabel.exe",
        "obabel",
    ]
    for c in candidates:
        try:
            r = subprocess.run([c, "--version"], capture_output=True, timeout=5)
            if r.returncode == 0 or b"Open Babel" in (r.stdout + r.stderr):
                return c
        except Exception as e:
            log.debug("_find_obabel: candidato %r indisponivel: %s", c, e)
    return "obabel"  # fallback to PATH


def ensure_default_compound():
    """
    Ensure the 2-SCN-NQ compound exists in the database.
    Call at startup to guarantee the base compound is always present.
    """
    name = "2-SCN-NQ"
    cdir = get_ligand_dir(name)
    meta_path = cdir / "metadata.json"

    if not meta_path.exists():
        cdir.mkdir(parents=True, exist_ok=True)
        with open(meta_path, "w", encoding='utf-8') as f:
            json.dump(DEFAULT_COMPOUND, f, indent=2, ensure_ascii=False)

    # Ensure PDBQT exists
    pdbqt = cdir / "2-SCN-NQ.pdbqt"
    if not pdbqt.exists():
        ok, msg = generate_pdbqt(DEFAULT_COMPOUND["smiles"], pdbqt, cdir / "2-SCN-NQ.mol2")
        return ok, msg
    return True, f"[OK] 2-SCN-NQ ready"


def populate_library(overwrite: bool = False) -> list[dict]:
    """
    Popula o banco com a biblioteca built-in de naftoquinonas antiparasitárias
    (NQ_LIBRARY). Cria metadata.json para cada composto; não sobrescreve
    entradas existentes a menos que overwrite=True. Não gera 3D/PDBQT
    (prepare_3d=False) para ser rápido — rodar prep separadamente se necessário.
    Retorna lista de {name, status} para exibir no dashboard.
    """
    results = []
    for entry in NQ_LIBRARY:
        name  = entry["name"]
        cdir  = get_ligand_dir(name)
        mpath = cdir / "metadata.json"
        if mpath.exists() and not overwrite:
            results.append({"name": name, "status": "já existe — pulado"})
            continue
        try:
            meta = add_compound(
                name=name,
                smiles=entry["smiles"],
                display_name=entry.get("display_name", name),
                notes=entry.get("notes", ""),
                prepare_3d=False,   # rápido; gerar 3D separadamente se necessário
            )
            results.append({"name": name, "status": "adicionado"})
        except Exception as e:
            results.append({"name": name, "status": f"ERRO: {e}"})
    return results


def clean_orphan_dirs() -> list[str]:
    """Remove subpastas em LIGAND_DIR sem metadata.json válido (ex: _object_Object_)."""
    removed = []
    if not LIGAND_DIR.exists():
        return removed
    for d in list(LIGAND_DIR.iterdir()):
        if not d.is_dir():
            continue
        mpath = d / "metadata.json"
        if not mpath.exists():
            try:
                import shutil
                shutil.rmtree(d)
                removed.append(d.name)
            except Exception as e:
                log.warning("clean_orphan_dirs: falha ao remover %s: %s", d, e)
    return removed


# =============================================================================
# SESSION MANAGEMENT
# =============================================================================

def make_session_tag() -> str:
    """Return a timestamp string for session subdirectories: YYYY-MM-DD_HH-MM-SS"""
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def get_session_dirs(base_dir: str | Path, session_tag: str, compound_name: str) -> dict:
    """
    Return a dict of timestamped output directories for a pipeline run.
    All per-session outputs live under figures/{compound}/{session}/ — a single
    canonical location.  The old docking/{compound}/{session}/ sub-tree was
    removed (it was never populated in practice and duplicated the layout).

    Structure:
      <base>/
        figures/          {compound}/{session}/   ← ALL per-session outputs
        manuscript_texts/ {compound}/{session}/
        logs/             {compound}/{session}/
        results/          {compound}/             ← cumulative, not per-session
        proteins/                                 ← shared receptors
        ligand/           {compound}/             ← prepared ligand files
    """
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    sess   = session_tag

    dirs = {
        "ligand":   base / "ligand"   / safe_c,
        "proteins": base / "proteins",
        "results":  base / "results"  / safe_c,
        "figures":  base / "figures"  / safe_c / sess,
        "logs":     base / "logs"     / safe_c / sess,
        "texts":    base / "manuscript_texts" / safe_c / sess,
        # "docking" alias kept for backward-compat callers; points to figures/
        "docking":  base / "figures"  / safe_c / sess,
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def list_sessions(base_dir: str | Path, compound_name: str) -> list[dict]:
    """
    List all past sessions for a compound, most recent first.
    Returns list of {tag, compound, figures_dir, docking_dir, texts_dir,
                     logs_dir, results_file}.
    Only looks in figures/ (single canonical location).
    """
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    sessions = []

    parent = base / "figures" / safe_c
    if parent.exists():
        for sess_dir in sorted(parent.iterdir(), reverse=True):
            if sess_dir.is_dir() and len(sess_dir.name) == 19:
                tag = sess_dir.name
                figures_dir = str(sess_dir)
                sessions.append({
                    "tag":          tag,
                    "compound":     compound_name,
                    "figures_dir":  figures_dir,
                    "docking_dir":  figures_dir,   # alias — same location
                    "texts_dir":    str(base / "manuscript_texts" / safe_c / tag),
                    "logs_dir":     str(base / "logs" / safe_c / tag),
                    "results_file": str(base / "results" / safe_c / "results.json"),
                })
    return sorted(sessions, key=lambda x: x["tag"], reverse=True)


def get_latest_results(base_dir: str | Path, compound_name: str) -> dict | None:
    """Return the most recent results.json for a compound, or None."""
    base   = Path(base_dir)
    safe_c = _safe_name(compound_name)
    rf     = base / "results" / safe_c / "results.json"
    if rf.exists():
        with open(rf) as f:
            return json.load(f)
    return None


# =============================================================================
# HELPERS
# =============================================================================

def _count_heavy_atoms(smiles: str) -> int:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol else 0
    except Exception:
        return 0


def _count_rotatable_bonds(smiles: str) -> int:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        return rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else 0
    except Exception:
        return 0


def _generate_svg(smiles: str, svg_path: Path) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem import Draw
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        mol = Draw.PrepareMolForDrawing(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_path.write_text(drawer.GetDrawingText())
        return True
    except Exception:
        return False


# =============================================================================
# CLI INTERFACE
# =============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compound Database Manager")
    sub = parser.add_subparsers(dest="cmd")

    # list
    sub.add_parser("list", help="List all compounds")

    # add
    p_add = sub.add_parser("add", help="Add a new compound")
    p_add.add_argument("--name",    required=True, help="Short ID (e.g. MY-COMPOUND-1)")
    p_add.add_argument("--smiles",  required=True)
    p_add.add_argument("--formula", default="")
    p_add.add_argument("--display", default="")
    p_add.add_argument("--notes",   default="")

    # info
    p_info = sub.add_parser("info", help="Show compound info")
    p_info.add_argument("name")

    # sessions
    p_sess = sub.add_parser("sessions", help="List past sessions for a compound")
    p_sess.add_argument("name")

    # init
    sub.add_parser("init", help="Ensure 2-SCN-NQ default compound exists")

    args = parser.parse_args()

    if args.cmd == "list":
        compounds = list_compounds()
        if not compounds:
            print("No compounds found. Run: python compound_db.py init")
        else:
            print(f"\n{'Name':<20} {'Formula':<16} {'MW':>8}  {'SMILES[:40]'}")
            print("-"*80)
            for c in compounds:
                smiles_preview = (c.get('smiles','') or '')[:40]
                pdbqt_ok = "[PDBQT]" if c.get("_has_pdbqt") else "[no PDBQT]"
                print(f"{c['name']:<20} {c.get('formula','?'):<16} {c.get('mw',0):>8.2f}  {smiles_preview} {pdbqt_ok}")

    elif args.cmd == "add":
        print(f"Adding compound: {args.name}")
        meta = add_compound(args.name, args.smiles, formula=args.formula,
                            display_name=args.display, notes=args.notes)
        print(f"  Directory: {meta['_dir']}")
        print(f"  Formula:   {meta['formula']}")
        print(f"  MW:        {meta['mw']}")
        print(f"  PDBQT:     {'OK' if meta['_has_pdbqt'] else 'MISSING'}")
        for msg in meta.get("_status_msgs", []):
            print(f"  {msg}")

    elif args.cmd == "info":
        c = load_compound(args.name)
        if c:
            for k, v in c.items():
                if not k.startswith("_"):
                    print(f"  {k:<20}: {v}")
        else:
            print(f"Compound '{args.name}' not found.")

    elif args.cmd == "sessions":
        sessions = list_sessions("docking_project", args.name)
        if not sessions:
            print(f"No sessions found for {args.name}")
        else:
            for s in sessions:
                print(f"  {s['tag']}  →  {s['figures_dir']}")

    elif args.cmd == "init":
        ok, msg = ensure_default_compound()
        print(msg)

    else:
        parser.print_help()
