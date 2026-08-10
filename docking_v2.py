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

# Module-level logger (used by module-level functions such as _convert_to_pdbqt).
# run_full_pipeline may rebind a local `log` to its own logger; this default
# ensures module functions never hit NameError when called outside that scope.
log = logging.getLogger("docking_v2")

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
# FONTE DE VERDADE: registry/targets/*.json (ver registry.py e docs/REGISTRY.md).
#
# O dicionario abaixo deixou de ser a fonte e passou a ser REDE DE SEGURANCA:
# so entra em uso se o registry estiver ausente ou corrompido, e nesse caso a
# plataforma avisa alto em vez de seguir em silencio. Manter as duas fontes em
# divergencia foi exatamente o que produziu tres PDBs errados no fallback de
# alvos do dashboard; o teste tests/test_registry_parity.py existe para impedir
# que isso se repita.
#
# Removido na Fase 9, quando o registry tiver rodagem suficiente.
_TARGET_DB_FALLBACK = {
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


def load_target_db(campaign_id=None):
    """Carrega os alvos do registry, caindo no fallback embutido se preciso.

    A campanha default vem de registry/program.json, nao de um literal aqui --
    e o que permite a mesma funcao servir naftoquinona em tripanotiona redutase
    e droga do citoesqueleto em actina sem tocar no codigo.
    """
    try:
        import registry

        campaign = campaign_id or registry.load_program().get("default_campaign", "SUB01")
        db = registry.as_legacy_target_db(campaign)
        if not db:
            raise registry.RegistryError(f"campanha '{campaign}' nao declara nenhum alvo")
        return db
    except Exception as exc:  # noqa: BLE001 - queremos degradar, mas ruidosamente
        import warnings

        warnings.warn(
            f"[docking_v2] registry indisponivel ({exc.__class__.__name__}: {exc}); "
            "usando o fallback embutido de alvos. Os resultados continuam validos para o "
            "SUB01, mas nenhuma campanha nova sera visivel ate o registry voltar.",
            RuntimeWarning,
            stacklevel=2,
        )
        return dict(_TARGET_DB_FALLBACK)


TARGET_DB = load_target_db()


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
                    if len(atoms) == 1:  # single atom = potential metal ion
                        metal_name = atoms[0].element.strip()
                        # Only count as metal when the real element is an
                        # essential metal; otherwise ions like CL/NA/K/BR
                        # (short resnames) polluted the QC metals report.
                        if (metal_name.upper() in ESSENTIAL_METALS
                                and metal_name not in report["metals_found"]):
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


def ensure_pdb_format(structure_path, output_dir, name_hint="structure"):
    """Garante um arquivo .pdb a partir de qualquer formato de estrutura aceito
    pelo registry (hoje: .pdb passa direto; .cif/.mmcif e convertido).

    clean_protein() usa Bio.PDB.PDBParser (formato PDB, colunas fixas de
    texto) -- parsear um .cif com ele NAO levanta excecao, mas produz uma
    Structure com ZERO atomos: falha silenciosa, confirmado por teste manual
    antes desta funcao existir (T-0021, Fase 0.3). Alvos com fonte AlphaFold/
    local (ex.: SUB04) sao .cif; alvos com PDB ID historico (SUB01) ja sao
    .pdb e retornam sem tocar em nada.
    """
    structure_path = Path(structure_path)
    suffix = structure_path.suffix.lower()
    if suffix == ".pdb":
        return structure_path
    if suffix in (".cif", ".mmcif"):
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(name_hint, str(structure_path))
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name_hint}_from_cif.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_path))
        return out_path
    raise ValueError(f"formato de estrutura desconhecido: {suffix!r} ({structure_path})")


def clean_protein(pdb_path, target_info, output_dir):
    """
    Clean protein keeping essential cofactors and metals.
    Returns: (clean_pdb_path, pdbqt_path, grid_center, qc_report)
    """
    from Bio.PDB import PDBParser, PDBIO, Select

    # Havia aqui um lookup reverso por PDB para descobrir o id do alvo. Era
    # codigo morto (o resultado nunca era usado) e uma bomba-relogio: estoura
    # IndexError com alvo sem PDB -- modelo AlphaFold, que e a regra fora do
    # painel original -- e devolve o alvo errado se dois alvos compartilharem
    # PDB (mesmo cristal, dois sitios). O id agora vem do proprio target_info.
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

    # Nome do arquivo de saida: PDB ID quando existe; alvos alphafold/local
    # (structure.pdb_id=None, ex.: os 5 do SUB04) usam o id do alvo -- sem
    # isso, processar 2+ alvos AlphaFold no mesmo output_dir colide em
    # "None_clean.pdb"/"None_receptor.pdbqt", sobrescrevendo em silencio.
    name_key = target_info.get("pdb") or target_info.get("id") or "unknown"

    # Save clean PDB
    clean_path = Path(output_dir) / f"{name_key}_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(clean_path), SmartSelect())

    # Convert to PDBQT
    pdbqt_path = Path(output_dir) / f"{name_key}_receptor.pdbqt"
    pdbqt_ok, pdbqt_msg = _convert_to_pdbqt(clean_path, pdbqt_path)

    return clean_path, pdbqt_path, grid_center, pdbqt_ok, pdbqt_msg


def _calculate_grid_center(structure, target_info):
    """Calculate grid center based on strategy"""
    from Bio.PDB import Selection
    import numpy as np

    strategy = target_info.get("grid_strategy", "ligand_center")

    # Strategy 0: centro EXPLICITO, ja resolvido em outro lugar (nao computado
    # aqui). Usado quando o alvo nao tem ligante co-cristalizado proprio (ex.:
    # modelos AlphaFold puros, caso do SUB04) -- o centro veio de uma
    # superposicao estrutural contra um complexo de referencia, feita fora
    # desta funcao (ver T-0021 Fase D) e gravada no registry como dado, nao
    # recomputada aqui. ERRO CLARO se a estrategia pede isto mas o numero nao
    # veio -- NUNCA cai no ultimo recurso (centro de massa da proteina
    # inteira), que seria um valor plausivel e silenciosamente errado.
    if strategy == "explicit_center":
        center = target_info.get("box_center_A")
        if not center:
            raise ValueError(
                "grid_strategy='explicit_center' mas target_info sem "
                "'box_center_A' -- o site precisa declarar o centro (nunca "
                "cair no centro de massa da proteina por engano)"
            )
        return tuple(round(float(c), 2) for c in center)

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
                        # Exact residue-identity match (resname+num), NOT substring:
                        # substring made "CYS52" also match "CYS520"/"CYS523", pulling
                        # wrong CA atoms into the grid-center average.
                        if any(fb.replace(" ", "") == full_name for fb in fallback_res):
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

    # Motivos de falha por rota. None = rota nao tentada; string = motivo real
    # (excecao ou saida invalida). Assim a mensagem final distingue "nao tentado"
    # de "tentado mas produziu saida pequena/invalida".
    _err_ob = None
    _err_rdk = None

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
        _err_ob = "saida <=500 bytes ou arquivo ausente"
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
        except Exception as e:
            log.debug("prepare receptor: SanitizeMol parcial falhou (%s); seguindo", e)

        # Calcular cargas Gasteiger
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except Exception as e:
            log.debug("prepare receptor: ComputeGasteigerCharges falhou (%s)", e)

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
        except Exception as e:
            log.debug("prepare receptor: deteccao de aromaticos falhou (%s)", e)

        # Acceptor N atoms (NA): aromatic N in heterocycle, ou N sem H ligado
        acceptor_n = set()
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N":
                num_h = atom.GetTotalNumHs()
                if num_h == 0 or atom.GetIsAromatic():
                    acceptor_n.add(atom.GetIdx())

        # Ler PDB original para preservar nomes de atomos / residuos / chains.
        # IMPORTANTE: NAO assumir alinhamento posicional 1:1 entre a lista de
        # atomos do RDKit e as linhas ATOM/HETATM do PDB. RDKit COLAPSA atomos
        # com indicador de conformacao alternativa (altLoc) — ex.: fica so com a
        # conformacao 'A' — resultando em MENOS atomos que o numero de linhas.
        # Como os atomos descartados ficam no MEIO do arquivo, o pareamento por
        # indice (pdb_lines[i] <-> rdkit_atoms[i]) desloca cargas/tipos para o
        # atomo ERRADO a partir do primeiro altLoc (medido: ~58% dos atomos
        # corrompidos em 1D3H; mismatch confirmado em 1D3H/2XU3/3C61/3DJJ/3KKU).
        # Pareamos por IDENTIDADE via numero de serie do atomo (cols 7-11), que
        # o RDKit preserva em GetPDBResidueInfo().GetSerialNumber().
        pdb_by_serial = {}
        pdb_line_count = 0
        with open(str(pdb_path), "r") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    pdb_line_count += 1
                    try:
                        serial = int(line[6:11])
                    except ValueError:
                        continue
                    pdb_by_serial.setdefault(serial, line)

        # Construir PDBQT receptor
        # Formato Vina:
        # ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67    -0.066 N
        rdkit_atoms = list(mol.GetAtoms())
        out_lines = ["REMARK  Generated by docking_v2.py (RDKit fallback)\n"]
        matched = 0
        unmatched = 0

        for atom in rdkit_atoms:
            idx = atom.GetIdx()

            # Localizar a linha PDB original por IDENTIDADE (numero de serie).
            # Sem correspondencia confiavel, NUNCA escrever uma linha (evita
            # atribuir carga/tipo ao atomo errado).
            info = atom.GetPDBResidueInfo()
            pdb_line = None
            if info is not None:
                pdb_line = pdb_by_serial.get(info.GetSerialNumber())
            if pdb_line is None:
                unmatched += 1
                continue
            matched += 1

            # rstrip("\r\n") cobre arquivos PDB com CRLF (Windows); apenas "\n"
            # deixaria \r residual e quebraria o slice [:66] mais abaixo.
            pdb_line = pdb_line.rstrip("\r\n")
            # Padronizar comprimento da linha (PDB tem 80 cols)
            if len(pdb_line) < 80:
                pdb_line = pdb_line + " " * (80 - len(pdb_line))

            # Pegar carga do RDKit
            try:
                q = float(atom.GetDoubleProp("_GasteigerCharge"))
                if q != q:  # NaN
                    q = 0.0
            except Exception:
                q = 0.0

            # Tipo AutoDock
            elem = atom.GetSymbol().upper()
            at = elem_to_at.get(elem, elem[:2])
            # Refinar tipo
            if elem == "C" and idx in aromatic_atoms:
                at = "A"
            elif elem == "N" and idx in acceptor_n:
                at = "NA"
            elif elem == "H":
                # H polar (HD) so se ligado a N/O/S
                neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
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

        # Guarda de seguranca: se o pareamento por identidade nao resolveu a
        # grande maioria dos atomos do RDKit, a saida seria um receptor vazio ou
        # gravemente incompleto. Abortar a rota (levanta -> capturado abaixo,
        # aciona a proxima rota / mensagem de erro) em vez de gravar um PDBQT
        # silenciosamente truncado/desalinhado.
        if matched == 0 or matched < 0.5 * len(rdkit_atoms):
            raise ValueError(
                f"Rota B RDKit: pareamento por identidade (serial) insuficiente "
                f"- {matched} de {len(rdkit_atoms)} atomos RDKit pareados "
                f"({unmatched} sem serial correspondente; {pdb_line_count} linhas "
                f"ATOM/HETATM). Abortando para nao gravar receptor desalinhado."
            )

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
            return True, f"RDKit fallback OK: {matched} atomos (pareados por serial)"
        _err_rdk = "saida <=500 bytes ou arquivo ausente"
    except Exception as e_rdk:
        _err_rdk = str(e_rdk)

    # Se ambas falharam:
    return False, (f"OpenBabel: {_err_ob if _err_ob is not None else 'nao tentado'} | "
                   f"RDKit: {_err_rdk if _err_rdk is not None else 'nao tentado'}")


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

def prepare_ligand_pdbqt(smiles, compound_name, ligand_root_dir, log=None):
    """Prepara o PDBQT de um ligante a partir do SMILES, com cache e cascata
    de rotas (OpenBabel -> RDKit+Meeko). Extraida de run_full_pipeline() (a
    Fase 4 previa isto -- 'run_full_pipeline ganha campaign_id/ligand_id' --
    mas o corpo continuava com SMILES/nome do SUB01 hardcoded ate agora), para
    ser reusavel por qualquer campanha/ligante via VinaEngine.prepare().

    Busca em 2 localizacoes (padrao do projeto) antes de gerar:
      1) <ligand_root_dir>/<compound_name>/<compound_name>.pdbqt (subpasta -- padrao compound_db)
      2) <ligand_root_dir>/<compound_name>.pdbqt                  (flat -- legado)

    Retorna o Path do PDBQT em sucesso, None em falha (log ja registra o motivo
    -- nunca levanta, para o chamador decidir como abortar).
    """
    log = log or logging.getLogger("pipeline")
    ligand_root_dir = Path(ligand_root_dir)

    candidate_paths = [
        ligand_root_dir / compound_name / f"{compound_name}.pdbqt",
        ligand_root_dir / f"{compound_name}.pdbqt",
    ]
    for cand in candidate_paths:
        if cand.exists() and cand.stat().st_size > 100:
            log.info(f"  [OK] Ligand exists: {cand}")
            return cand

    # Nao encontrou -- gerar. Usa subpasta (padrao do projeto)
    ligand_dir = ligand_root_dir / compound_name
    ligand_dir.mkdir(parents=True, exist_ok=True)
    ligand_pdbqt = ligand_dir / f"{compound_name}.pdbqt"

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
        mol.write("mol2", str(ligand_dir / f"{compound_name}.mol2"), overwrite=True)
        log.info(f"  [OK] Ligand prepared via OpenBabel: {ligand_pdbqt}")
        return ligand_pdbqt
    except Exception as e_ob:
        # Rota B: RDKit + Meeko -- tentada em QUALQUER falha da Rota A, nao so
        # ImportError (pacote ausente). Achado real (T-0021 Fase E): OpenBabel
        # pode estar INSTALADO mas mal configurado no ambiente (aqui:
        # "smi is not a recognised Open Babel format", tipico de
        # BABEL_DATADIR/BABEL_LIBDIR ausente) -- isso levantava uma excecao
        # generica que a versao anterior desta funcao NAO tratava como "tentar
        # a Rota B", so registrava o erro e desistia.
        log.info(f"  [INFO] Rota A (OpenBabel) falhou ({e_ob!r}) — tentando Rota B (RDKit + Meeko)")
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
            mol2_path = ligand_dir / f"{compound_name}.mol2"
            Chem.MolToMolFile(rdmol, str(mol2_path.with_suffix(".mol")))
            log.info(f"  [OK] Ligand prepared via RDKit+Meeko: {ligand_pdbqt}")
            return ligand_pdbqt
        except Exception as e_rdk:
            log.error(f"  [ERRO] Rota B (RDKit+Meeko) falhou: {e_rdk}")
            log.error(f"  [ERRO] Para gerar o ligante, instale uma dessas opcoes:")
            log.error(f"         conda install -c conda-forge openbabel  (Rota A)")
            log.error(f"         pip install meeko                       (Rota B - precisa rdkit)")
            log.error(f"         OU mantenha o env 'datascience' ativo (ja tem tudo):")
            log.error(f"         conda activate datascience")
            return None


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

    ligand_pdbqt = prepare_ligand_pdbqt(smiles, compound_safe, dirs["ligand"], log=log)
    if ligand_pdbqt is None:
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
        except Exception as e:
            log.debug("receptor QC: checagem de PDBQT corrompido falhou (%s)", e)

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

        # Save results JSON.
        # Grava os DOIS nomes durante a transicao: 'results.json' e a convencao
        # canonica nova, 'docking_results.json' e o que leitores antigos esperam.
        # A escrita dupla cai na Fase 9, quando nao houver mais leitor antigo.
        import json
        results_json = dirs["results"] / "docking_results.json"
        try:
            results_json.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
            log.info(f"\n  [OK] Results saved: {results_json}")
        except Exception as e:
            log.warning(f"  [AVISO] Could not save JSON: {e}")

        canonical = dirs["results"] / "results.json"
        try:
            canonical.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
            log.info(f"  [OK] Results saved (canonico): {canonical}")
        except Exception as e:
            log.warning(f"  [AVISO] Could not save canonical JSON: {e}")

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
