#!/usr/bin/env python3
"""
=============================================================================
PROTEIN PREP + QC MODULE - Docking Pipeline v2
Corrected based on literature:
- FAD MUST be retained in TR structures
- NADPH/NADP+ MUST be retained in PTR1
- Metal cofactors (Zn, Fe, Mn) retained
- Grid calculated from co-crystallized ligand center
- L. amazonensis TR has no crystal structure (use 2JK6 + homology note)
=============================================================================
"""
import os
import sys
import json
import math
import shutil
import subprocess
import logging
import requests
import argparse
from pathlib import Path
from datetime import datetime

# =====================================================================
# COFACTORS AND ESSENTIAL HETEROATOMS TO KEEP
# =====================================================================
# These residue names are KEPT during protein cleaning
ESSENTIAL_COFACTORS = {
    "FAD",   # Flavin adenine dinucleotide - essential for TR, DHODH
    "FMN",   # Flavin mononucleotide - essential for DHODH
    "NAD",   # Nicotinamide adenine dinucleotide
    "NAP",   # NADPH
    "NDP",   # NADP+
    "NAI",   # NAD+ (alternate code)
    "FDA",   # FAD (alternate code)
}

ESSENTIAL_METALS = {
    "ZN",    # Zinc
    "FE",    # Iron (Fe-SOD)
    "FE2",   # Iron (Fe2+)
    "MN",    # Manganese (MnSOD)
    "MG",    # Magnesium
    "CA",    # Calcium
    "CO",    # Cobalt
    "CU",    # Copper
}

# Water and common buffer molecules to REMOVE
REMOVE_RESIDUES = {
    "HOH", "WAT", "DOD",   # Water
    "SO4", "SUL",           # Sulfate
    "PO4", "PHO",           # Phosphate
    "GOL",                  # Glycerol
    "EDO",                  # Ethylene glycol
    "ACT",                  # Acetate
    "CIT",                  # Citrate
    "PEG",                  # PEG
    "DMS", "DMF",           # Solvents
    "BME",                  # Beta-mercaptoethanol
    "TRS",                  # Tris buffer
    "CL",                   # Chloride ion
    "IOD",                  # Iodide
    "MPD",                  # MPD
}

# =====================================================================
# TARGET DATABASE v2 - Literature-informed
# =====================================================================
TARGET_DB = {
    # --- T. cruzi ---
    "TcTR": {
        "pdb": "1BZL",
        "organism": "Trypanosoma cruzi",
        "target": "Trypanothione Reductase",
        "category": "parasite",
        "keep_cofactors": ["FAD"],
        "ligand_resnames": ["TS2", "TS4"],  # Trypanothione
        "grid_strategy": "ligand_center",    # center on co-crystal ligand
        "grid_fallback_residues": ["CYS52", "CYS57", "HIS461"],
        "grid_size": (25, 25, 25),
        "human_ortholog": "hGR",
        "notes": "FAD is ESSENTIAL. Grid centered on trypanothione binding site.",
        "literature": "Krauth-Siegel 2003; PDB has FAD + TS2 co-crystal.",
    },
    "TcCruzain": {
        "pdb": "3KKU",
        "organism": "Trypanosoma cruzi",
        "target": "Cruzain",
        "category": "parasite",
        "keep_cofactors": [],
        "ligand_resnames": ["K11", "VNS"],  # vinylsulfone inhibitors
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": ["CYS25", "HIS159", "ASN175"],
        "grid_size": (22, 22, 22),
        "human_ortholog": "hCatL",
        "notes": "No cofactors needed. Grid on catalytic triad.",
        "literature": "McKerrow 2009; co-crystal with K11.",
    },
    "TcFeSOD": {
        "pdb": "4DVH",
        "organism": "Trypanosoma cruzi",
        "target": "Fe-Superoxide Dismutase",
        "category": "parasite",
        "keep_cofactors": ["FE", "FE2"],
        "ligand_resnames": [],  # apo structure
        "grid_strategy": "metal_center",     # center on Fe atom
        "grid_fallback_residues": ["HIS26", "HIS73", "ASP156", "HIS160"],
        "grid_size": (22, 22, 22),
        "human_ortholog": "hMnSOD",
        "notes": "Fe is ESSENTIAL. Apo structure - grid on Fe active site.",
        "literature": "Mateo 2008.",
    },
    # --- Leishmania ---
    "LiTR": {
        "pdb": "2JK6",
        "organism": "Leishmania infantum (proxy for L. amazonensis, 95% identity)",
        "target": "Trypanothione Reductase",
        "category": "parasite",
        "keep_cofactors": ["FAD"],
        "ligand_resnames": ["TS2", "TS4", "FAD"],
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": ["CYS52", "CYS57", "HIS461"],
        "grid_size": (25, 25, 25),
        "human_ortholog": "hGR",
        "notes": (
            "NO crystal structure exists for L. amazonensis TR. "
            "2JK6 (L. infantum) is the standard proxy (95% identity). "
            "Camargo et al. (Sci Rep, 2023) used this approach for L. amazonensis. "
            "FAD is ESSENTIAL - grid centered on flavin ring near active site."
        ),
        "literature": "Baiocco 2009; Camargo/Sci Rep 2023 for L. amazonensis.",
    },
    "LmPTR1": {
        "pdb": "2QHX",
        "organism": "Leishmania major",
        "target": "Pteridine Reductase 1 (PTR1)",
        "category": "parasite",
        "keep_cofactors": ["NAP", "NDP", "FAD"],  # NADPH essential
        "ligand_resnames": ["6QA", "MTX", "TAQ", "DHB"],
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": ["PHE113", "TYR194", "SER111"],
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "notes": (
            "NADPH is ESSENTIAL - it forms part of the active site. "
            "2QHX has co-crystal with pteridine inhibitor (best for grid). "
            "Alternative: 1E7W (methotrexate complex)."
        ),
        "literature": "Gourley 2001; Dawson 2006; Hardy 2005.",
    },
    "LmDHODH": {
        "pdb": "3C61",
        "organism": "Leishmania donovani",
        "target": "Dihydroorotate Dehydrogenase",
        "category": "parasite",
        "keep_cofactors": ["FMN", "FAD"],
        "ligand_resnames": ["ORO", "OA"],  # orotate
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": [],
        "grid_size": (25, 25, 25),
        "human_ortholog": "hDHODH",
        "notes": "FMN essential. Naphthoquinones may compete as ubiquinone analogs.",
        "literature": "Pinheiro 2013.",
    },
    # --- Human controls ---
    "hGR": {
        "pdb": "3DJJ",
        "organism": "Homo sapiens",
        "target": "Glutathione Reductase",
        "category": "human",
        "keep_cofactors": ["FAD"],
        "ligand_resnames": ["GSH", "GTX", "GDS"],  # glutathione
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": ["CYS58", "CYS63"],
        "grid_size": (25, 25, 25),
        "human_ortholog": None,
        "notes": "FAD essential. Counterpart of TR.",
        "literature": "Karplus & Schulz 1987.",
    },
    "hCatL": {
        "pdb": "2XU3",
        "organism": "Homo sapiens",
        "target": "Cathepsin L",
        "category": "human",
        "keep_cofactors": [],
        "ligand_resnames": [],
        "grid_strategy": "residue_center",
        "grid_fallback_residues": ["CYS25", "HIS163", "ASN187"],
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "notes": "Counterpart of cruzain.",
        "literature": "Turk 2012.",
    },
    "hDHODH": {
        "pdb": "1D3H",
        "organism": "Homo sapiens",
        "target": "DHODH",
        "category": "human",
        "keep_cofactors": ["FMN"],
        "ligand_resnames": ["BRE"],  # brequinar
        "grid_strategy": "ligand_center",
        "grid_fallback_residues": [],
        "grid_size": (25, 25, 25),
        "human_ortholog": None,
        "notes": "FMN essential.",
        "literature": "Liu 2000.",
    },
    "hMnSOD": {
        "pdb": "1LUV",
        "organism": "Homo sapiens",
        "target": "Mn-Superoxide Dismutase",
        "category": "human",
        "keep_cofactors": ["MN"],
        "ligand_resnames": [],
        "grid_strategy": "metal_center",
        "grid_fallback_residues": ["HIS26", "HIS74", "ASP159", "HIS163"],
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "notes": "Mn essential. Counterpart of Fe-SOD.",
        "literature": "Borgstahl 1992.",
    },
}


# =====================================================================
# QC + PROTEIN PREPARATION
# =====================================================================

def download_pdb(pdb_id, output_dir):
    """Download PDB from RCSB"""
    path = Path(output_dir) / f"{pdb_id}.pdb"
    if path.exists():
        return path, "[OK] Already exists"

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        path.write_text(r.text)
        return path, "[OK] Downloaded"
    except Exception as e:
        return None, f"[ERRO] Download failed: {e}"


def analyze_pdb(pdb_path, target_info):
    """Quality control: analyze PDB contents"""
    from Bio.PDB import PDBParser

    report = {
        "file": str(pdb_path),
        "pdb_id": target_info["pdb"],
        "chains": [],
        "residue_count": 0,
        "water_count": 0,
        "cofactors_found": [],
        "metals_found": [],
        "ligands_found": [],
        "other_het": [],
        "expected_cofactors": target_info.get("keep_cofactors", []),
        "expected_ligands": target_info.get("ligand_resnames", []),
        "warnings": [],
        "status": "OK",
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(pdb_path))

    chains = set()
    for model in structure:
        for chain in model:
            chains.add(chain.id)
            for residue in chain:
                het_flag = residue.id[0]
                resname = residue.get_resname().strip()

                if het_flag == " ":
                    report["residue_count"] += 1
                elif het_flag == "W":
                    report["water_count"] += 1
                elif resname in ESSENTIAL_COFACTORS:
                    if resname not in report["cofactors_found"]:
                        report["cofactors_found"].append(resname)
                elif resname in ESSENTIAL_METALS or len(resname) <= 2:
                    # Metals are typically 1-2 char residue names
                    atoms = list(residue.get_atoms())
                    if len(atoms) == 1:  # single atom = metal ion
                        metal_name = atoms[0].element.strip()
                        if metal_name and metal_name not in report["metals_found"]:
                            report["metals_found"].append(metal_name)
                elif resname in REMOVE_RESIDUES:
                    pass  # buffer, will be removed
                else:
                    if resname not in report["other_het"]:
                        report["other_het"].append(resname)
                    if resname in target_info.get("ligand_resnames", []):
                        if resname not in report["ligands_found"]:
                            report["ligands_found"].append(resname)

    report["chains"] = sorted(chains)

    # Warnings
    expected_cof = target_info.get("keep_cofactors", [])
    for cof in expected_cof:
        if cof not in report["cofactors_found"] and cof not in report["metals_found"]:
            if cof in ESSENTIAL_METALS:
                if cof not in [m.upper() for m in report["metals_found"]]:
                    report["warnings"].append(
                        f"Expected cofactor {cof} NOT FOUND in PDB!")
            else:
                report["warnings"].append(
                    f"Expected cofactor {cof} NOT FOUND in PDB!")

    if not report["ligands_found"] and target_info.get("ligand_resnames"):
        report["warnings"].append(
            "No expected co-crystallized ligand found. "
            "Grid will use fallback residue-based centering.")

    if report["warnings"]:
        report["status"] = "WARNING"

    return report


def clean_protein(pdb_path, target_info, output_dir):
    """
    Clean protein keeping essential cofactors and metals.
    Returns: (clean_pdb_path, pdbqt_path, grid_center, qc_report)
    """
    from Bio.PDB import PDBParser, PDBIO, Select

    target_id = [k for k, v in TARGET_DB.items() if v["pdb"] == target_info["pdb"]][0]
    keep_cofactors = set(target_info.get("keep_cofactors", []))
    keep_metals = ESSENTIAL_METALS

    class SmartSelect(Select):
        """Keep protein + cofactors + metals, remove water + buffer"""
        def accept_residue(self, residue):
            het_flag = residue.id[0]
            resname = residue.get_resname().strip()

            # Standard amino acids - always keep
            if het_flag == " ":
                return True

            # Water - always remove
            if het_flag == "W" or resname in REMOVE_RESIDUES:
                return False

            # Essential cofactors for this target
            if resname in keep_cofactors:
                return True

            # Metal ions - keep if single atom
            atoms = list(residue.get_atoms())
            if len(atoms) == 1:
                elem = atoms[0].element.strip().upper()
                if elem in keep_metals or resname.upper() in keep_metals:
                    return True

            # Co-crystallized ligand - REMOVE (we dock our own compound)
            if resname in target_info.get("ligand_resnames", []):
                return False

            # Unknown HETATM - remove by default
            return False

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(pdb_path))

    # Calculate grid center BEFORE cleaning (need ligand coordinates)
    grid_center = _calculate_grid_center(structure, target_info)

    # Save clean PDB
    clean_path = Path(output_dir) / f"{target_info['pdb']}_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(clean_path), SmartSelect())

    # Convert to PDBQT
    pdbqt_path = Path(output_dir) / f"{target_info['pdb']}_receptor.pdbqt"
    pdbqt_ok, pdbqt_msg = _convert_to_pdbqt(clean_path, pdbqt_path)

    return clean_path, pdbqt_path, grid_center, pdbqt_ok, pdbqt_msg


def _calculate_grid_center(structure, target_info):
    """Calculate grid center based on strategy"""
    from Bio.PDB import Selection
    import numpy as np

    strategy = target_info.get("grid_strategy", "ligand_center")

    # Strategy 1: Center on co-crystallized ligand
    if strategy == "ligand_center":
        ligand_resnames = target_info.get("ligand_resnames", [])
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in ligand_resnames and resname not in ESSENTIAL_COFACTORS:
                        atoms = list(residue.get_atoms())
                        if atoms:
                            coords = np.array([a.get_vector().get_array() for a in atoms])
                            center = coords.mean(axis=0)
                            return tuple(round(float(c), 2) for c in center)

    # Strategy 2: Center on metal cofactor
    if strategy == "metal_center":
        for model in structure:
            for chain in model:
                for residue in chain:
                    atoms = list(residue.get_atoms())
                    if len(atoms) == 1:
                        elem = atoms[0].element.strip().upper()
                        expected = [c.upper() for c in target_info.get("keep_cofactors", [])]
                        if elem in expected or elem in ESSENTIAL_METALS:
                            coord = atoms[0].get_vector().get_array()
                            return tuple(round(float(c), 2) for c in coord)

    # Strategy 3 / Fallback: Center on named residues
    fallback_res = target_info.get("grid_fallback_residues", [])
    if fallback_res:
        all_coords = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":
                        resname = residue.get_resname().strip()
                        resnum = str(residue.id[1])
                        full_name = f"{resname}{resnum}"
                        if any(fb.replace(" ", "") in full_name for fb in fallback_res):
                            for atom in residue:
                                if atom.name == "CA":
                                    all_coords.append(atom.get_vector().get_array())
        if all_coords:
            coords = np.array(all_coords)
            center = coords.mean(axis=0)
            return tuple(round(float(c), 2) for c in center)

    # Last resort: center of mass of protein
    all_coords = []
    for atom in structure.get_atoms():
        all_coords.append(atom.get_vector().get_array())
    if all_coords:
        coords = np.array(all_coords)
        center = coords.mean(axis=0)
        return tuple(round(float(c), 2) for c in center)

    return (0.0, 0.0, 0.0)


def _convert_to_pdbqt(pdb_path, pdbqt_path):
    """Convert PDB to PDBQT receptor format.

    Tenta multiplas rotas em ordem de preferencia:
      Rota A: OpenBabel (precisa BABEL_DATADIR configurado)
      Rota B: Conversao manual PDB->PDBQT (atribui Gasteiger via RDKit)
      Rota C: Meeko (se disponivel)

    Vina receptors devem ter apenas ATOM/HETATM/REMARK/END (sem ROOT/BRANCH).
    """

    # ── Rota A: OpenBabel ──────────────────────────────────────────────────
    try:
        from openbabel import openbabel as ob
        from openbabel import pybel

        # Verificar se o formato PDB esta disponivel
        formats = pybel.informats
        if "pdb" not in formats:
            raise ImportError(
                f"OpenBabel sem formato PDB ({len(formats)} formatos disponiveis). "
                f"Provavel BABEL_DATADIR nao configurado."
            )

        mol = next(pybel.readfile("pdb", str(pdb_path)))
        mol.addh()

        charge_model = ob.OBChargeModel.FindType("gasteiger")
        if charge_model:
            charge_model.ComputeCharges(mol.OBMol)

        temp_path = str(pdbqt_path) + ".tmp"
        mol.write("pdbqt", temp_path, overwrite=True)

        # Strip ligand tags do PDBQT
        bad_tags = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
        kept_lines = 0
        removed_lines = 0
        with open(temp_path, "r") as fin, open(str(pdbqt_path), "w") as fout:
            for line in fin:
                if line.strip().startswith(bad_tags):
                    removed_lines += 1
                    continue
                fout.write(line)
                kept_lines += 1
        Path(temp_path).unlink(missing_ok=True)

        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 500:
            return True, f"OpenBabel OK: {kept_lines} linhas, {removed_lines} tags removidas"
    except Exception as e_ob:
        _err_ob = str(e_ob)

    # ── Rota B: PDB->PDBQT manual via RDKit + tabela AutoDock ────────────────
    # Esta rota nao precisa de OpenBabel datafiles; usa RDKit para cargas
    # Gasteiger e mapeia tipos AutoDock baseado em elemento + hibridizacao.
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Ler PDB com RDKit (sem sanitizar, para tolerar PDBs nao-canonicos)
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
        if mol is None:
            raise ValueError("RDKit MolFromPDBFile retornou None")
        try:
            Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES))
        except Exception:
            pass

        # Calcular cargas Gasteiger
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except Exception:
            pass

        # Mapeamento simples de elemento -> tipo AutoDock para receptor
        # (proteinas tem elementos limitados: C, N, O, S, H, P, e ocasionalmente metais)
        elem_to_at = {
            "C":  "C",   # carbono alifatico
            "N":  "N",   # nitrogenio (sera ajustado para NA se aceitador)
            "O":  "OA",  # oxigenio (geralmente aceitador em proteinas)
            "S":  "SA",  # enxofre (Cys, Met)
            "P":  "P",   # fosforo
            "H":  "HD",  # hidrogenio polar (apenas se ligado a heteroatomo)
            "F":  "F",
            "CL": "Cl",
            "BR": "Br",
            "I":  "I",
            "FE": "Fe",
            "MN": "Mn",
            "MG": "Mg",
            "ZN": "Zn",
            "CA": "Ca",
            "K":  "K",
            "NA": "Na",
        }

        # Conjunto de aromaticos (ajusta C->A, N->NA, etc.)
        aromatic_atoms = set()
        try:
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    aromatic_atoms.add(atom.GetIdx())
        except Exception:
            pass

        # Acceptor N atoms (NA): aromatic N in heterocycle, ou N sem H ligado
        acceptor_n = set()
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N":
                num_h = atom.GetTotalNumHs()
                if num_h == 0 or atom.GetIsAromatic():
                    acceptor_n.add(atom.GetIdx())

        # Ler PDB original para preservar nomes de atomos / residuos / chains
        pdb_lines = []
        with open(str(pdb_path), "r") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    pdb_lines.append(line)

        # Construir PDBQT receptor
        # Formato Vina:
        # ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67    -0.066 N
        rdkit_atoms = list(mol.GetAtoms())
        out_lines = ["REMARK  Generated by docking_v2.py (RDKit fallback)\n"]
        n_atoms = min(len(pdb_lines), len(rdkit_atoms))

        for i in range(n_atoms):
            # rstrip("\r\n") cobre arquivos PDB com CRLF (Windows); apenas "\n"
            # deixaria \r residual e quebraria o slice [:66] mais abaixo.
            pdb_line = pdb_lines[i].rstrip("\r\n")
            # Padronizar comprimento da linha (PDB tem 80 cols)
            if len(pdb_line) < 80:
                pdb_line = pdb_line + " " * (80 - len(pdb_line))

            # Pegar carga do RDKit
            try:
                q = float(rdkit_atoms[i].GetDoubleProp("_GasteigerCharge"))
                if q != q:  # NaN
                    q = 0.0
            except Exception:
                q = 0.0

            # Tipo AutoDock
            elem = rdkit_atoms[i].GetSymbol().upper()
            at = elem_to_at.get(elem, elem[:2])
            # Refinar tipo
            if elem == "C" and i in aromatic_atoms:
                at = "A"
            elif elem == "N" and i in acceptor_n:
                at = "NA"
            elif elem == "H":
                # H polar (HD) so se ligado a N/O/S
                neighbors = [n.GetSymbol() for n in rdkit_atoms[i].GetNeighbors()]
                if any(n in ("N", "O", "S") for n in neighbors):
                    at = "HD"
                else:
                    continue  # pular H apolar (Vina nao precisa)

            # Formato PDBQT Vina (79 chars):
            #   cols 1-66 = ATOM/HETATM record do PDB original (com occ/tempf)
            #   cols 67-70 = 4 espacos
            #   cols 71-76 = charge formato %+6.3f (sinal explicito)
            #   col 77 = espaco
            #   cols 78-79 = atom type AutoDock (left-aligned, 2 chars)
            # Exemplo: "ATOM   1  O  UNL 1       0.900  -2.725   0.027  0.00  0.00    +0.000 OA"
            charge_str = f"{q:+6.3f}"  # ex: "+0.000", "-0.320"
            line_out = pdb_line[:66] + "    " + charge_str + " " + f"{at:<2}\n"
            out_lines.append(line_out)

        out_lines.append("END\n")

        # Escrever
        with open(str(pdbqt_path), "w") as fout:
            fout.writelines(out_lines)

        # ── Self-test do formato PDBQT ────────────────────────────────────
        # Vina rejeita silenciosamente charges fora das cols 71-76. Validar
        # imediatamente para falhar rapido em vez de deixar a docking explodir
        # alvo a alvo. Bug historico: `q:>8.3f` ou rstrip('\n') com CRLF
        # deslocava o campo charge — agora travamos cedo.
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 500:
            try:
                with open(str(pdbqt_path), "r") as _fin:
                    for _ln in _fin:
                        if _ln.startswith(("ATOM", "HETATM")):
                            _ch = _ln[70:76]
                            float(_ch)  # ValueError se mal-formatado
                            break
            except (ValueError, IndexError):
                return False, (
                    f"Self-test PDBQT falhou: campo charge cols 71-76 "
                    f"corrompido. Investigue _convert_to_pdbqt rota B."
                )
            return True, f"RDKit fallback OK: {n_atoms} atomos"
    except Exception as e_rdk:
        _err_rdk = str(e_rdk)

    # Se ambas falharam:
    return False, (f"OpenBabel: {_err_ob if '_err_ob' in dir() else 'nao tentado'} | "
                   f"RDKit: {_err_rdk if '_err_rdk' in dir() else 'nao tentado'}")


# =====================================================================
# DOCKING EXECUTION
# =====================================================================

def find_vina():
    """Find vina executable"""
    # Check PATH
    vina = shutil.which("vina")
    if vina:
        return vina
    # Check common locations
    for candidate in [
        Path("vina.exe"), Path("vina"),
        Path("docking_project/vina.exe"), Path("docking_project/vina"),
        Path("./tools/vina.exe"), Path("./tools/vina"),
    ]:
        if candidate.exists():
            return str(candidate.resolve())
    return None


def run_docking(receptor_pdbqt, ligand_pdbqt, output_pdbqt,
                grid_center, grid_size, exhaustiveness=32, vina_path=None):
    """Run AutoDock Vina docking"""
    vina_cmd = vina_path or find_vina()
    if not vina_cmd:
        return None, "Vina executable not found"

    cx, cy, cz = grid_center
    sx, sy, sz = grid_size

    cmd = [
        vina_cmd,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "10",
        "--energy_range", "3",
        "--out", str(output_pdbqt),
    ]

    try:
        # Log the full command for debugging
        cmd_str = " ".join(cmd)
        print(f"    [DEBUG] Command: {cmd_str[:200]}...")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)

        # Capture both stdout and stderr
        full_output = result.stdout
        if result.stderr:
            full_output += "\n--- STDERR ---\n" + result.stderr

        if result.returncode != 0:
            return None, f"Vina exit code {result.returncode}.\nSTDOUT: {result.stdout[:200]}\nSTDERR: {result.stderr[:200]}"

        # Parse energies from stdout
        affinities = []
        reading = False
        for line in result.stdout.splitlines():
            if "-----+------------" in line:
                reading = True
                continue
            if reading and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinities.append(float(parts[1]))
                    except ValueError:
                        break

        return affinities, full_output
    except subprocess.TimeoutExpired:
        return None, "Vina timed out (>20 min)"
    except Exception as e:
        return None, str(e)


# =====================================================================
# PYMOL LAUNCHER
# =====================================================================

def find_pymol():
    """Find PyMOL executable"""
    import platform
    candidates = []
    if platform.system() == "Windows":
        candidates = [
            Path(os.environ.get("PROGRAMFILES", "")) / "PyMOL" / "PyMOLWin.exe",
            Path(os.environ.get("PROGRAMFILES(X86)", "")) / "PyMOL" / "PyMOLWin.exe",
            Path(os.environ.get("LOCALAPPDATA", "")) / "Schrodinger" / "PyMOL2" / "PyMOLWin.exe",
            Path("C:/PyMOL/PyMOLWin.exe"),
        ]
    else:
        candidates = [Path("/usr/bin/pymol"), Path("/usr/local/bin/pymol")]

    # Also check PATH
    pymol_path = shutil.which("pymol")
    if pymol_path:
        return pymol_path

    for c in candidates:
        if c.exists():
            return str(c)

    return None


def open_in_pymol(files, pymol_path=None):
    """Open files in PyMOL"""
    pymol = pymol_path or find_pymol()
    if not pymol:
        return False, "PyMOL not found"

    cmd = [pymol] + [str(f) for f in files]
    try:
        subprocess.Popen(cmd)
        return True, "[OK] PyMOL opened"
    except Exception as e:
        return False, f"[ERRO] {e}"


# =====================================================================
# FULL PIPELINE
# =====================================================================

def run_full_pipeline(selected_targets, base_dir="docking_project",
                     exhaustiveness=32, open_pymol=False):
    """Run complete pipeline with QC"""

    base = Path(base_dir)
    dirs = {k: base / k for k in ["ligand", "proteins", "docking", "results", "figures", "logs"]}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Logging
    import io
    log_file = dirs["logs"] / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    stream = (io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
              if hasattr(sys.stdout, 'buffer') else sys.stdout)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file, encoding='utf-8'),
                  logging.StreamHandler(stream)]
    )
    log = logging.getLogger("pipeline")

    log.info("=" * 60)
    log.info("DOCKING PIPELINE v2 - with QC + cofactor retention")
    log.info(f"Targets: {', '.join(selected_targets)}")
    log.info("=" * 60)

    # --- LIGAND PREP ---
    log.info("\n--- STEP 1: Ligand Preparation ---")
    smiles = "O=C1C=C(SC#N)C(=O)c2ccccc21"
    compound_safe = "2-SCN-NQ"  # Nome do composto base do estudo

    # Buscar ligante PDBQT em multiplas localizacoes (padrao do projeto):
    #  1) docking_project/ligand/2-SCN-NQ/2-SCN-NQ.pdbqt   (subpasta - padrao compound_db)
    #  2) docking_project/ligand/2-SCN-NQ.pdbqt             (flat - legado)
    candidate_paths = [
        dirs["ligand"] / compound_safe / f"{compound_safe}.pdbqt",
        dirs["ligand"] / f"{compound_safe}.pdbqt",
    ]
    ligand_pdbqt = None
    for cand in candidate_paths:
        if cand.exists() and cand.stat().st_size > 100:
            ligand_pdbqt = cand
            log.info(f"  [OK] Ligand exists: {ligand_pdbqt}")
            break

    if ligand_pdbqt is None:
        # Nao encontrou — tentar gerar. Usa subpasta (padrao do projeto)
        ligand_dir = dirs["ligand"] / compound_safe
        ligand_dir.mkdir(parents=True, exist_ok=True)
        ligand_pdbqt = ligand_dir / f"{compound_safe}.pdbqt"

        # Rota A: OpenBabel (preferido, gera PDBQT diretamente com cargas Gasteiger)
        try:
            from openbabel import pybel
            mol = pybel.readstring("smi", smiles)
            mol.addh()
            mol.make3D(forcefield="mmff94", steps=500)
            ff = pybel._forcefields["mmff94"]
            if ff.Setup(mol.OBMol):
                ff.ConjugateGradients(1000, 1.0e-6)
                ff.GetCoordinates(mol.OBMol)
            mol.write("pdbqt", str(ligand_pdbqt), overwrite=True)
            mol.write("mol2", str(ligand_dir / f"{compound_safe}.mol2"), overwrite=True)
            log.info(f"  [OK] Ligand prepared via OpenBabel: {ligand_pdbqt}")
        except ImportError:
            # Rota B: RDKit + Meeko (fallback sem openbabel)
            log.info("  [INFO] OpenBabel ausente — tentando Rota B (RDKit + Meeko)")
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                from meeko import MoleculePreparation, PDBQTWriterLegacy

                rdmol = Chem.MolFromSmiles(smiles)
                if rdmol is None:
                    raise ValueError("SMILES invalido")
                rdmol = Chem.AddHs(rdmol)
                ret = AllChem.EmbedMolecule(rdmol, AllChem.ETKDGv3())
                if ret == -1:
                    ret = AllChem.EmbedMolecule(rdmol)
                if ret == -1:
                    raise ValueError("EmbedMolecule falhou")
                AllChem.MMFFOptimizeMolecule(rdmol)

                prep = MoleculePreparation()
                prep.prepare(rdmol)
                pdbqt_text = PDBQTWriterLegacy.write_string(prep.setup)[0]
                ligand_pdbqt.write_text(pdbqt_text, encoding="utf-8")

                # Tambem salvar mol2 via RDKit
                mol2_path = ligand_dir / f"{compound_safe}.mol2"
                Chem.MolToMolFile(rdmol, str(mol2_path.with_suffix(".mol")))
                log.info(f"  [OK] Ligand prepared via RDKit+Meeko: {ligand_pdbqt}")
            except Exception as e_rdk:
                log.error(f"  [ERRO] Rota B (RDKit+Meeko) falhou: {e_rdk}")
                log.error(f"  [ERRO] Para gerar o ligante, instale uma dessas opcoes:")
                log.error(f"         conda install -c conda-forge openbabel  (Rota A)")
                log.error(f"         pip install meeko                       (Rota B - precisa rdkit)")
                log.error(f"         OU mantenha o env 'datascience' ativo (ja tem tudo):")
                log.error(f"         conda activate datascience")
                return False
        except Exception as e:
            log.error(f"  [ERRO] OpenBabel falhou: {e}")
            return False

    # --- PROTEIN PREP + QC ---
    log.info("\n--- STEP 2: Protein Preparation + Quality Control ---")

    # Check for corrupt PDBQT files from old pipeline (ROOT tag issue)
    for existing_pdbqt in dirs["proteins"].glob("*_receptor.pdbqt"):
        try:
            first_lines = existing_pdbqt.read_text(errors='ignore')[:2000]
            if "ROOT" in first_lines or "BRANCH" in first_lines:
                log.warning(f"  [AVISO] Corrupt PDBQT detected: {existing_pdbqt.name}")
                log.warning(f"    Contains ROOT/BRANCH tags (ligand format, not receptor)")
                log.warning(f"    Deleting to re-generate correctly...")
                existing_pdbqt.unlink()
        except Exception:
            pass

    all_results = {}
    vina_path = find_vina()

    for tid in selected_targets:
        if tid not in TARGET_DB:
            log.warning(f"  [AVISO] Unknown target: {tid}")
            continue

        info = TARGET_DB[tid]
        pdb_id = info["pdb"]
        log.info(f"\n  === {tid}: {info['target']} ({info['organism']}) ===")

        # Download
        pdb_path, msg = download_pdb(pdb_id, dirs["proteins"])
        log.info(f"    Download: {msg}")
        if not pdb_path:
            continue

        # QC Analysis
        log.info(f"    Running quality control...")
        qc = analyze_pdb(pdb_path, info)
        log.info(f"    Chains: {qc['chains']}")
        log.info(f"    Residues: {qc['residue_count']}")
        log.info(f"    Waters: {qc['water_count']} (will be removed)")
        log.info(f"    Cofactors found: {qc['cofactors_found']} (KEEPING)")
        log.info(f"    Metals found: {qc['metals_found']} (KEEPING)")
        log.info(f"    Ligands found: {qc['ligands_found']}")
        log.info(f"    Other HET: {qc['other_het']} (removing)")
        if qc["warnings"]:
            for w in qc["warnings"]:
                log.warning(f"    [AVISO] {w}")
        log.info(f"    QC Status: {qc['status']}")

        # Clean + grid calculation
        log.info(f"    Cleaning protein (keeping cofactors)...")
        clean_pdb, pdbqt, grid_center, pdbqt_ok, pdbqt_msg = clean_protein(
            pdb_path, info, dirs["proteins"])

        if pdbqt_ok:
            log.info(f"    [OK] Clean PDB: {clean_pdb}")
            log.info(f"    [OK] PDBQT: {pdbqt} ({pdbqt_msg})")
        else:
            log.warning(f"    [AVISO] PDBQT conversion failed: {pdbqt_msg}")
            log.warning(f"    Prepare manually in Chimera: DockPrep -> Save as PDBQT")

        log.info(f"    Grid center: {grid_center}")
        log.info(f"    Grid size: {info['grid_size']}")

        # Notes
        if info.get("notes"):
            log.info(f"    Note: {info['notes']}")

        # --- DOCKING ---
        log.info(f"    Running docking...")

        # Verify grid center is valid
        if grid_center == (0.0, 0.0, 0.0):
            log.error(f"    [ERRO] Grid center is (0,0,0) - could not determine active site!")
            log.error(f"    Check if PDB has co-crystallized ligand for grid centering.")
            continue

        if not vina_path:
            log.error("    [ERRO] Vina not found!")
            log.error("    Download from: https://github.com/ccsb-scripps/AutoDock-Vina/releases")
            continue

        receptor = pdbqt if pdbqt_ok else None
        if not receptor or not receptor.exists():
            log.error(f"    [ERRO] No receptor PDBQT for {tid}")
            log.error(f"    Expected: {pdbqt}")
            continue

        # Verify receptor file is valid
        rec_size = receptor.stat().st_size
        if rec_size < 500:
            log.error(f"    [ERRO] Receptor PDBQT too small ({rec_size} bytes)")
            log.error(f"    Delete {receptor} and re-run protein preparation.")
            continue

        # Verify ligand exists and is valid
        if not ligand_pdbqt.exists() or ligand_pdbqt.stat().st_size < 100:
            log.error(f"    [ERRO] Ligand PDBQT missing or invalid")
            continue

        output = dirs["docking"] / f"result_{tid}.pdbqt"
        log.info(f"    Receptor: {receptor} ({rec_size:,} bytes)")
        log.info(f"    Ligand: {ligand_pdbqt} ({ligand_pdbqt.stat().st_size:,} bytes)")
        log.info(f"    Grid: center={grid_center} size={info['grid_size']}")

        affinities, vina_output = run_docking(
            receptor, ligand_pdbqt, output,
            grid_center, info["grid_size"],
            exhaustiveness, vina_path
        )

        # Save vina log
        log_dock = dirs["docking"] / f"log_{tid}.txt"
        if vina_output:
            log_dock.write_text(vina_output)

        if affinities:
            best_dG = affinities[0]
            ki = math.exp(best_dG / (1.987e-3 * 298.15)) * 1e6
            log.info(f"    [OK] Best dG = {best_dG:.1f} kcal/mol | Ki = {ki:.3f} uM")
            log.info(f"    [OK] All modes: {[f'{a:.1f}' for a in affinities[:5]]}")

            # Verify output file
            if output.exists() and output.stat().st_size > 100:
                log.info(f"    [OK] Output: {output} ({output.stat().st_size:,} bytes)")
            else:
                log.warning(f"    [AVISO] Output file missing or empty: {output}")

            all_results[tid] = {
                "pdb": pdb_id,
                "organism": info["organism"],
                "target": info["target"],
                "category": info["category"],
                "best_dG": best_dG,
                "ki_uM": ki,
                "all_modes": affinities[:5],
                "grid_center": list(grid_center),
                "grid_size": list(info["grid_size"]),
                "cofactors_kept": qc["cofactors_found"],
                "qc_status": qc["status"],
            }
        else:
            err_msg = vina_output[:500] if vina_output else "no output"
            log.error(f"    [ERRO] Docking failed for {tid}!")
            log.error(f"    Vina output: {err_msg}")

    # --- RESULTS ---
    log.info("\n--- STEP 3: Results ---")

    n_attempted = len(selected_targets)
    n_success = len(all_results)
    n_failed = n_attempted - n_success

    if all_results:
        log.info(f"\n  Docking succeeded for {n_success}/{n_attempted} targets:")
        log.info(f"\n  {'Target':<12} {'Organism':<35} {'dG':<10} {'Ki (uM)':<10}")
        log.info(f"  {'-'*12} {'-'*35} {'-'*10} {'-'*10}")
        for tid, data in sorted(all_results.items(), key=lambda x: x[1]["best_dG"]):
            log.info(f"  {tid:<12} {data['organism']:<35} "
                     f"{data['best_dG']:<10.1f} {data['ki_uM']:<10.3f}")

        # Save results JSON
        import json
        results_json = dirs["results"] / "docking_results.json"
        try:
            results_json.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
            log.info(f"\n  [OK] Results saved: {results_json}")
        except Exception as e:
            log.warning(f"  [AVISO] Could not save JSON: {e}")

    if n_failed > 0:
        log.warning(f"\n  {n_failed} target(s) failed - see log above for details.")

    if not all_results:
        log.error(f"\n  [ERRO] NO DOCKING RESULTS GENERATED!")
        log.error(f"  All {n_attempted} targets failed.")
        log.error(f"  Common causes:")
        log.error(f"    1. PDBQT receptor files are empty or corrupt")
        log.error(f"    2. Grid center is at (0,0,0) - wrong position")
        log.error(f"    3. Vina.exe not found or wrong version")
        log.error(f"    4. Ligand PDBQT is invalid")
        log.error(f"  Run: python diagnostic.py  for detailed checks")

    # --- PYMOL ---
    if open_pymol:
        pymol = find_pymol()
        if pymol:
            log.info(f"\n  Opening PyMOL...")
            files = [dirs["proteins"] / f"{TARGET_DB[t]['pdb']}_clean.pdb"
                     for t in selected_targets if t in TARGET_DB]
            open_in_pymol(files, pymol)
        else:
            log.info("  [AVISO] PyMOL not found")

    log.info("\n" + "=" * 60)
    log.info("PIPELINE COMPLETE")
    log.info("=" * 60)
    return True


# =====================================================================
# CLI
# =====================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Docking Pipeline v2 - QC + Cofactors")
    parser.add_argument("--targets", type=str, default=None,
                       help="Comma-separated target IDs (e.g., TcTR,LiTR,hGR)")
    parser.add_argument("--all", action="store_true", help="All targets")
    parser.add_argument("--parasite", action="store_true", help="All parasite targets")
    parser.add_argument("--exhaustiveness", type=int, default=32)
    parser.add_argument("--pymol", action="store_true", help="Open results in PyMOL")
    parser.add_argument("--list", action="store_true", help="List targets")
    parser.add_argument("--qc-only", action="store_true",
                       help="Only run QC analysis, no docking")

    args = parser.parse_args()

    if args.list:
        print("\nAvailable targets:")
        print(f"  {'ID':<12} {'PDB':<6} {'Organism':<40} {'Target'}")
        print(f"  {'-'*12} {'-'*6} {'-'*40} {'-'*30}")
        for tid, info in TARGET_DB.items():
            print(f"  {tid:<12} {info['pdb']:<6} {info['organism']:<40} {info['target']}")
        sys.exit(0)

    if args.all:
        targets = list(TARGET_DB.keys())
    elif args.parasite:
        targets = [k for k, v in TARGET_DB.items() if v["category"] == "parasite"]
    elif args.targets:
        targets = [t.strip() for t in args.targets.split(",")]
    else:
        print("Use --targets TcTR,LiTR,hGR  or  --all  or  --list")
        sys.exit(1)

    run_full_pipeline(targets, exhaustiveness=args.exhaustiveness,
                     open_pymol=args.pymol)
