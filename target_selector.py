#!/usr/bin/env python3
"""
=============================================================================
SELETOR INTERATIVO DE ALVOS PARA DOCKING MOLECULAR
Composto 2: 2-Tiocianato-1,4-Naftoquinona (2-SCN-NQ)

Uso:
    python target_selector.py              # Menu interativo
    python target_selector.py --targets 1,3,5,10   # Selecao direta
    python target_selector.py --all-parasite        # Todos os parasitarios
    python target_selector.py --all                 # Todos (parasito + humano)
=============================================================================
"""

import os
import sys
import json
import math
import shutil
import subprocess
import argparse
import logging
import requests
from pathlib import Path
from datetime import datetime

# =============================================================================
# BANCO DE DADOS DE ALVOS
# =============================================================================

TARGET_DATABASE = {

    # =========================================================================
    # TRYPANOSOMA CRUZI
    # =========================================================================
    1: {
        "id": "TcTR",
        "pdb": "1BZL",
        "resolution": "2.40 A",
        "organism": "T. cruzi",
        "target_name": "Tripanotiona Redutase (TcTR)",
        "target_name_en": "Trypanothione Reductase",
        "category": "parasite",
        "justification": (
            "Alvo classico e validado. Ausente em mamiferos (humanos usam "
            "glutationa redutase). Diretamente coerente com dados de estresse "
            "oxidativo - inibicao da TcTR compromete defesa antioxidante "
            "tiol-dependente do parasito, explicando aumento de ROS."
        ),
        "cofactor": "FAD, NADPH",
        "ligand_in_pdb": "Tripanotiona (TS2)",
        "grid_center": (15.0, -2.0, -18.0),
        "grid_size": (25, 25, 25),
        "human_ortholog": "hGR",
        "references": "Krauth-Siegel et al., 2003, Antioxid Redox Signal",
    },
    2: {
        "id": "TcCruzain",
        "pdb": "3KKU",
        "resolution": "1.28 A",
        "organism": "T. cruzi",
        "target_name": "Cruzaina (cisteino-protease)",
        "target_name_en": "Cruzain (cysteine protease)",
        "category": "parasite",
        "justification": (
            "Principal cisteino-protease de T. cruzi. Essencial para "
            "sobrevivencia, nutricao e evasao imune do parasito. Co-cristal "
            "com inibidor covalente disponivel (excelente para validacao)."
        ),
        "cofactor": "Nenhum",
        "ligand_in_pdb": "Inibidor vinilsulfona (K11)",
        "grid_center": (5.0, 10.0, 15.0),
        "grid_size": (22, 22, 22),
        "human_ortholog": "hCatL",
        "references": "McKerrow et al., 2009, Bioorg Med Chem",
    },
    3: {
        "id": "TcFeSOD",
        "pdb": "4DVH",
        "resolution": "1.50 A",
        "organism": "T. cruzi",
        "target_name": "Ferro-Superoxido Dismutase (Fe-SOD)",
        "target_name_en": "Iron Superoxide Dismutase (Fe-SOD)",
        "category": "parasite",
        "justification": (
            "Diretamente relacionada ao estresse oxidativo. Fe-SOD e a "
            "principal defesa contra O2.- no parasito. Ausente em humanos "
            "(que usam Cu/Zn-SOD e Mn-SOD). Inibicao potencializa dano "
            "oxidativo da naftoquinona - sinergia com dados de ROS."
        ),
        "cofactor": "Fe2+/Fe3+",
        "ligand_in_pdb": "Nenhum (apo)",
        "grid_center": (30.0, 45.0, 25.0),
        "grid_size": (22, 22, 22),
        "human_ortholog": "hMnSOD",
        "references": "Mateo et al., 2008, Mol Biochem Parasitol",
    },
    4: {
        "id": "TcNDH2",
        "pdb": "5JWA",
        "resolution": "2.60 A",
        "organism": "T. cruzi",
        "target_name": "NADH desidrogenase tipo II (NDH2)",
        "target_name_en": "Type II NADH Dehydrogenase (NDH2)",
        "category": "parasite",
        "justification": (
            "NADH-fumarato redutase mitocondrial / NADH desidrogenase tipo II. "
            "Ausente em mamiferos. Componente da cadeia respiratoria alternativa "
            "de tripanossomatideos. Coerente com dados de estresse oxidativo - "
            "naftoquinonas podem desviar eletrons gerando ROS mitocondrial. "
            "NOTA: PDB 5JWA e de Caldalkalibacillus (melhor estrutura "
            "disponivel); usar como modelo por homologia se necessario."
        ),
        "cofactor": "FAD, UQ (ubiquinona)",
        "ligand_in_pdb": "UQ1",
        "grid_center": (25.0, 20.0, 15.0),
        "grid_size": (25, 25, 25),
        "human_ortholog": None,
        "references": "Verner et al., 2011, Int J Parasitol",
    },

    # =========================================================================
    # LEISHMANIA
    # =========================================================================
    5: {
        "id": "LiTR",
        "pdb": "2JK6",
        "resolution": "2.95 A",
        "organism": "Leishmania infantum",
        "target_name": "Tripanotiona Redutase (LiTR)",
        "target_name_en": "Trypanothione Reductase",
        "category": "parasite",
        "justification": (
            "Ortologo da TcTR em Leishmania. Mesmo mecanismo - essencial "
            "para manutencao do sistema tiol-dependente. Permite comparacao "
            "direta T. cruzi vs Leishmania no manuscrito."
        ),
        "cofactor": "FAD, NADPH",
        "ligand_in_pdb": "Tripanotiona",
        "grid_center": (10.5, 3.0, -14.0),
        "grid_size": (25, 25, 25),
        "human_ortholog": "hGR",
        "references": "Baiocco et al., 2009, Structure",
    },
    6: {
        "id": "LmTR",
        "pdb": "2W0H",
        "resolution": "2.60 A",
        "organism": "Leishmania major",
        "target_name": "Tripanotiona Redutase (LmTR)",
        "target_name_en": "Trypanothione Reductase (L. major)",
        "category": "parasite",
        "justification": (
            "Alternativa ao LiTR com especie diferente (L. major). "
            "Util se quiser mostrar conservacao do binding entre especies."
        ),
        "cofactor": "FAD, NADPH",
        "ligand_in_pdb": "Inibidor",
        "grid_center": (12.0, 5.0, -12.0),
        "grid_size": (25, 25, 25),
        "human_ortholog": "hGR",
        "references": "Khan, 2007, J Mol Graph Model",
    },
    7: {
        "id": "LmPTR1",
        "pdb": "1E7W",
        "resolution": "2.20 A",
        "organism": "Leishmania major",
        "target_name": "Pteridina Redutase 1 (PTR1)",
        "target_name_en": "Pteridine Reductase 1 (PTR1)",
        "category": "parasite",
        "justification": (
            "Alvo validado com excelentes estruturas cristalograficas. "
            "Essencial para metabolismo de pteridinas/folatos. Bypass da "
            "via DHFR - resistencia a antifolatos. Varias co-estruturas "
            "com inibidores disponiveis para validacao do protocolo."
        ),
        "cofactor": "NADPH",
        "ligand_in_pdb": "Methotrexato (MTX)",
        "grid_center": (25.0, 30.0, 15.0),
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "references": "Gourley et al., 2001, Nat Struct Biol",
    },
    8: {
        "id": "LmDHODH",
        "pdb": "3C61",
        "resolution": "2.13 A",
        "organism": "Leishmania donovani",
        "target_name": "Diidroorotato Desidrogenase (DHODH)",
        "target_name_en": "Dihydroorotate Dehydrogenase (DHODH)",
        "category": "parasite",
        "justification": (
            "Enzima da via de biossintese de pirimidinas. Alvo emergente "
            "para tripanossomatideos. Usa ubiquinona como aceptor de "
            "eletrons - naftoquinonas podem competir como analogos. "
            "Estruturalmente distinta da DHODH humana."
        ),
        "cofactor": "FMN, Ubiquinona",
        "ligand_in_pdb": "Orotato (ORO)",
        "grid_center": (-34.34, 21.12, -2.42),
        "grid_size": (25, 25, 25),
        "human_ortholog": "hDHODH",
        "references": "Pinheiro et al., 2013, Biochemistry",
    },

    # =========================================================================
    # ORTOLOGOS HUMANOS (CONTROLES DE SELETIVIDADE)
    # =========================================================================
    9: {
        "id": "hGR",
        "pdb": "3DJJ",
        "resolution": "1.05 A",
        "organism": "H. sapiens",
        "target_name": "Glutationa Redutase (hGR)",
        "target_name_en": "Glutathione Reductase",
        "category": "human",
        "justification": (
            "Contraparte humana da tripanotiona redutase. Usa glutationa "
            "em vez de tripanotiona. Diferenca estrutural no sitio ativo "
            "(mais estreito) e a base molecular da seletividade."
        ),
        "cofactor": "FAD, NADPH",
        "ligand_in_pdb": "FAD/NADPH (sem GSSG; grid no par catalitico Cys58/Cys63)",
        "grid_center": (20.48, 14.94, 25.71),
        "grid_size": (25, 25, 25),
        "human_ortholog": None,
        "references": "Karplus & Schulz, 1987, J Mol Biol",
    },
    10: {
        "id": "hCatL",
        "pdb": "2XU3",
        "resolution": "1.40 A",
        "organism": "H. sapiens",
        "target_name": "Catepsina L humana",
        "target_name_en": "Cathepsin L",
        "category": "human",
        "justification": (
            "Contraparte humana da cruzaina. Comparacao direta de "
            "seletividade entre cisteino-proteases parasitaria vs hospedeiro."
        ),
        "cofactor": "Nenhum",
        "ligand_in_pdb": "Inibidor",
        "grid_center": (10.0, 15.0, 20.0),
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "references": "Turk et al., 2012, Biochim Biophys Acta",
    },
    11: {
        "id": "hDHODH",
        "pdb": "1D3H",
        "resolution": "1.60 A",
        "organism": "H. sapiens",
        "target_name": "DHODH humana",
        "target_name_en": "Dihydroorotate Dehydrogenase",
        "category": "human",
        "justification": (
            "Contraparte humana da DHODH de Leishmania. Diferencas "
            "no sitio de ligacao da ubiquinona permitem seletividade."
        ),
        "cofactor": "FMN, Ubiquinona",
        "ligand_in_pdb": "Brequinar",
        "grid_center": (25.0, 20.0, 15.0),
        "grid_size": (25, 25, 25),
        "human_ortholog": None,
        "references": "Liu et al., 2000, Structure",
    },
    12: {
        "id": "hMnSOD",
        "pdb": "1LUV",
        "resolution": "1.85 A",
        "organism": "H. sapiens",
        "target_name": "Mn-Superoxido Dismutase (MnSOD)",
        "target_name_en": "Manganese Superoxide Dismutase",
        "category": "human",
        "justification": (
            "Contraparte humana da Fe-SOD de T. cruzi. Usa Mn em vez "
            "de Fe como cofator - diferenca no sitio ativo permite "
            "seletividade. Comparacao Fe-SOD vs MnSOD e classica."
        ),
        "cofactor": "Mn3+/Mn2+",
        "ligand_in_pdb": "Nenhum (apo)",
        "grid_center": (28.0, 32.0, 20.0),
        "grid_size": (22, 22, 22),
        "human_ortholog": None,
        "references": "Borgstahl et al., 1992, Cell",
    },
}

# Pares de seletividade (parasito --> humano)
SELECTIVITY_PAIRS = {
    "TcTR": "hGR",
    "LiTR": "hGR",
    "LmTR": "hGR",
    "TcCruzain": "hCatL",
    "TcFeSOD": "hMnSOD",
    "LmDHODH": "hDHODH",
}

# Configuracao do composto
COMPOUND = {
    "name": "2-SCN-NQ",
    "smiles": "O=C1C=C(SC#N)C(=O)c2ccccc21",
    "formula": "C11H5NO2S",
    "mw": 215.23,
}


# =============================================================================
# INTERFACE DE SELECAO
# =============================================================================

def display_menu():
    """Exibe menu de selecao de alvos"""

    print("\n" + "=" * 80)
    print("  SELETOR DE ALVOS - Docking Molecular: 2-SCN-NQ")
    print("=" * 80)

    # T. cruzi
    print("\n  +-------------------------------------------------------------+")
    print("  |  TRYPANOSOMA CRUZI                                         |")
    print("  +-------------------------------------------------------------+")
    for num, data in TARGET_DATABASE.items():
        if data["organism"].startswith("T. cruzi"):
            pdb = data["pdb"]
            name = data["target_name"]
            res = data["resolution"]
            print(f"  |  [{num:>2}]  {name:<42} PDB: {pdb}  ({res}) |")
    print("  +-------------------------------------------------------------+")

    # Leishmania
    print("\n  +-------------------------------------------------------------+")
    print("  |  LEISHMANIA spp.                                           |")
    print("  +-------------------------------------------------------------+")
    for num, data in TARGET_DATABASE.items():
        if "Leishmania" in data["organism"]:
            pdb = data["pdb"]
            name = data["target_name"]
            res = data["resolution"]
            print(f"  |  [{num:>2}]  {name:<42} PDB: {pdb}  ({res}) |")
    print("  +-------------------------------------------------------------+")

    # Humanos
    print("\n  +-------------------------------------------------------------+")
    print("  |  ORTOLOGOS HUMANOS (controles de seletividade)             |")
    print("  +-------------------------------------------------------------+")
    for num, data in TARGET_DATABASE.items():
        if data["category"] == "human":
            pdb = data["pdb"]
            name = data["target_name"]
            res = data["resolution"]
            orth = [k for k, v in SELECTIVITY_PAIRS.items() if v == data["id"]]
            pair_str = f"< vs {', '.join(orth)}" if orth else ""
            print(f"  |  [{num:>2}]  {name:<42} PDB: {pdb}  {pair_str:<10} |")
    print("  +-------------------------------------------------------------+")

    # Atalhos
    print("\n  Atalhos:")
    print("    [A]  Todos os alvos (parasitarios + humanos)")
    print("    [P]  Todos os parasitarios (T. cruzi + Leishmania)")
    print("    [T]  So T. cruzi (alvos 1-4)")
    print("    [L]  So Leishmania (alvos 5-8)")
    print("    [H]  So humanos (alvos 9-12)")
    print("    [S]  Selecao recomendada para manuscrito (TcTR + LiTR + hGR)")
    print("    [Q]  Sair")


def get_user_selection():
    """Obtem selecao do usuario"""
    print("\n" + "-" * 60)
    print("  Digite os numeros separados por virgula (ex: 1,5,9)")
    print("  ou use um atalho (A/P/T/L/H/S/Q):")
    print("-" * 60)

    choice = input("\n  >>> ").strip().upper()

    if choice == 'Q':
        print("\n  Saindo...")
        sys.exit(0)
    elif choice == 'A':
        return list(TARGET_DATABASE.keys())
    elif choice == 'P':
        return [n for n, d in TARGET_DATABASE.items() if d["category"] == "parasite"]
    elif choice == 'T':
        return [n for n, d in TARGET_DATABASE.items() if d["organism"].startswith("T. cruzi")]
    elif choice == 'L':
        return [n for n, d in TARGET_DATABASE.items() if "Leishmania" in d["organism"]]
    elif choice == 'H':
        return [n for n, d in TARGET_DATABASE.items() if d["category"] == "human"]
    elif choice == 'S':
        # Selecao recomendada: TcTR + LiTR + hGR (par principal de seletividade)
        return [1, 5, 9]
    else:
        try:
            nums = [int(x.strip()) for x in choice.split(",")]
            valid = [n for n in nums if n in TARGET_DATABASE]
            invalid = [n for n in nums if n not in TARGET_DATABASE]
            if invalid:
                print(f"\n  [AVISO] Alvos invalidos ignorados: {invalid}")
            return valid
        except ValueError:
            print("\n  [ERRO] Entrada invalida. Use numeros separados por virgula.")
            return []


def display_selection_summary(selected_nums):
    """Mostra resumo da selecao e identifica controles humanos necessarios"""
    print("\n" + "=" * 80)
    print("  ALVOS SELECIONADOS")
    print("=" * 80)

    selected = {n: TARGET_DATABASE[n] for n in selected_nums}
    parasite_targets = {n: d for n, d in selected.items() if d["category"] == "parasite"}
    human_targets = {n: d for n, d in selected.items() if d["category"] == "human"}

    # Mostrar selecionados
    for num, data in selected.items():
        icon = "[P]" if data["category"] == "parasite" else "[H]"
        print(f"\n  {icon} [{num}] {data['target_name']} ({data['organism']})")
        print(f"      PDB: {data['pdb']} | Resolucao: {data['resolution']}")
        print(f"      Cofator: {data['cofactor']}")
        print(f"      Ligante no PDB: {data['ligand_in_pdb']}")
        # Justificativa com wrapping simples
        just = data['justification']
        words = just.split()
        line = "      --> "
        for w in words:
            if len(line) + len(w) > 75:
                print(line)
                line = "        " + w
            else:
                line += " " + w if line.strip() else w
        if line.strip():
            print(line)

    # Sugerir controles humanos se nao incluidos
    suggested = set()
    for num, data in parasite_targets.items():
        orth = data.get("human_ortholog")
        if orth:
            # Verificar se o ortologo humano ja esta na selecao
            human_ids = {d["id"] for d in human_targets.values()}
            if orth not in human_ids:
                # Encontrar o numero do ortologo
                for hnum, hdata in TARGET_DATABASE.items():
                    if hdata["id"] == orth:
                        suggested.add((hnum, hdata["target_name"], orth,
                                      data["target_name"]))

    if suggested:
        print("\n" + "-" * 60)
        print("  [AVISO] CONTROLES HUMANOS RECOMENDADOS (nao selecionados):")
        print("-" * 60)
        for hnum, hname, hid, pname in suggested:
            print(f"    [{hnum}] {hname} < contraparte de {pname}")

        print(f"\n  Adicionar controles humanos automaticamente? (S/N)")
        add = input("  >>> ").strip().upper()
        if add in ('S', 'Y', 'SIM', 'YES', ''):
            for hnum, _, _, _ in suggested:
                if hnum not in selected_nums:
                    selected_nums.append(hnum)
            print(f"  [OK] Controles adicionados. Total: {len(selected_nums)} alvos")

    return selected_nums


def display_target_detail(num):
    """Exibe detalhes completos de um alvo"""
    if num not in TARGET_DATABASE:
        print(f"  [ERRO] Alvo {num} nao encontrado.")
        return

    d = TARGET_DATABASE[num]
    print(f"\n{'=' * 60}")
    print(f"  DETALHE: [{num}] {d['target_name']}")
    print(f"{'=' * 60}")
    print(f"  ID interno:     {d['id']}")
    print(f"  PDB:            {d['pdb']}")
    print(f"  Resolucao:      {d['resolution']}")
    print(f"  Organismo:      {d['organism']}")
    print(f"  Categoria:      {d['category']}")
    print(f"  Cofator:        {d['cofactor']}")
    print(f"  Ligante no PDB: {d['ligand_in_pdb']}")
    print(f"  Grid center:    {d['grid_center']}")
    print(f"  Grid size:      {d['grid_size']}")
    if d.get('human_ortholog'):
        print(f"  Ortologo humano: {d['human_ortholog']}")
    print(f"  Referencia:     {d['references']}")
    print(f"\n  Justificativa:")
    print(f"  {d['justification']}")


# =============================================================================
# EXECUCAO DO DOCKING PARA ALVOS SELECIONADOS
# =============================================================================

def run_docking_for_selection(selected_nums, base_dir="docking_project"):
    """Executa pipeline de docking para os alvos selecionados"""

    base = Path(base_dir)
    ligand_dir = base / "ligand"
    protein_dir = base / "proteins"
    docking_dir = base / "docking"
    results_dir = base / "results"
    figures_dir = base / "figures"
    logs_dir = base / "logs"

    for d in [ligand_dir, protein_dir, docking_dir, results_dir, figures_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = logs_dir / f"docking_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    import io
    stream_handler = logging.StreamHandler(
        io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if hasattr(sys.stdout, 'buffer') else sys.stdout
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file, encoding='utf-8'), stream_handler]
    )
    log = logging.getLogger(__name__)

    log.info("+======================================================+")
    log.info("|  DOCKING MOLECULAR - 2-SCN-NQ                        |")
    log.info(f"|  Alvos selecionados: {len(selected_nums):<33}|")
    log.info(f"|  Data: {datetime.now().strftime('%Y-%m-%d %H:%M'):<46}|")
    log.info("+======================================================+")

    # ----- STEP 1: Preparar ligante -----
    log.info("\n" + "=" * 60)
    log.info("STEP 1: Preparacao do Ligante")
    log.info("=" * 60)

    ligand_pdbqt = ligand_dir / f"{COMPOUND['name']}.pdbqt"
    ligand_sdf = ligand_dir / f"{COMPOUND['name']}.sdf"

    if not ligand_pdbqt.exists():
        try:
            from openbabel import pybel
            mol = pybel.readstring("smi", COMPOUND["smiles"])
            mol.addh()
            mol.make3D(forcefield="mmff94", steps=500)

            # Otimizar
            ff = pybel._forcefields["mmff94"]
            if ff.Setup(mol.OBMol):
                ff.ConjugateGradients(1000, 1.0e-6)
                ff.GetCoordinates(mol.OBMol)

            mol.write("sdf", str(ligand_sdf), overwrite=True)
            mol.write("pdbqt", str(ligand_pdbqt), overwrite=True)
            mol.write("mol2", str(ligand_dir / f"{COMPOUND['name']}.mol2"), overwrite=True)
            log.info(f"  [OK] Ligante preparado: {ligand_pdbqt}")
            log.info(f"    Atomos pesados: {mol.OBMol.NumHvyAtoms()}")
            log.info(f"    Rotatable bonds: {mol.OBMol.NumRotors()}")
        except ImportError:
            log.error("  [ERRO] Open Babel nao encontrado!")
            log.error("    pip install openbabel-wheel")
            log.error("    OU: conda install -c conda-forge openbabel")
            return False
    else:
        log.info(f"  [OK] Ligante ja existe: {ligand_pdbqt}")

    # ----- STEP 2: Baixar e preparar proteinas -----
    log.info("\n" + "=" * 60)
    log.info("STEP 2: Preparacao das Proteinas")
    log.info("=" * 60)

    for num in selected_nums:
        target = TARGET_DATABASE[num]
        pdb_id = target["pdb"]
        log.info(f"\n  --- [{num}] {target['target_name']} ({target['organism']}) ---")

        pdb_path = protein_dir / f"{pdb_id}.pdb"
        clean_path = protein_dir / f"{pdb_id}_clean.pdb"
        pdbqt_path = protein_dir / f"{pdb_id}_receptor.pdbqt"

        # Baixar
        if not pdb_path.exists():
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            log.info(f"    Baixando {pdb_id}...")
            try:
                r = requests.get(url, timeout=30)
                r.raise_for_status()
                with open(pdb_path, "w") as f:
                    f.write(r.text)
                log.info(f"    [OK] Baixado: {pdb_path}")
            except Exception as e:
                log.error(f"    [ERRO] Erro: {e}")
                continue
        else:
            log.info(f"    [OK] PDB ja existe")

        # Limpar e preparar
        if not pdbqt_path.exists():
            try:
                from Bio.PDB import PDBParser, PDBIO, Select

                class ProteinOnly(Select):
                    def accept_residue(self, residue):
                        return residue.id[0] == " "

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("prot", str(pdb_path))
                io = PDBIO()
                io.set_structure(structure)
                io.save(str(clean_path), ProteinOnly())

                # Converter para PDBQT
                from openbabel import pybel
                pmol = next(pybel.readfile("pdb", str(clean_path)))
                pmol.addh()
                pmol.write("pdbqt", str(pdbqt_path), overwrite=True)
                log.info(f"    [OK] Receptor PDBQT: {pdbqt_path}")

            except Exception as e:
                log.warning(f"    [AVISO] Preparacao automatica falhou: {e}")
                log.info(f"    --> Prepare manualmente no Chimera e salve como:")
                log.info(f"      {pdbqt_path}")
        else:
            log.info(f"    [OK] Receptor ja preparado")

    # ----- STEP 3: Docking -----
    log.info("\n" + "=" * 60)
    log.info("STEP 3: Docking Molecular (AutoDock Vina)")
    log.info("=" * 60)

    all_results = {}
    use_python_vina = False
    vina_executable = None

    try:
        from vina import Vina
        use_python_vina = True
        log.info("  Engine: Python binding do Vina")
    except ImportError:
        # Buscar executavel: PATH, pasta do projeto, pasta atual
        vina_executable = shutil.which("vina")
        if not vina_executable:
            for candidate in [base / "vina.exe", base / "vina",
                             Path("vina.exe"), Path("vina")]:
                if candidate.exists():
                    vina_executable = str(candidate.resolve())
                    break
        if vina_executable:
            log.info(f"  Engine: Executavel vina ({vina_executable})")
        else:
            log.error("  [ERRO] AutoDock Vina nao encontrado!")
            log.error("    Baixe de: https://github.com/ccsb-scripps/AutoDock-Vina/releases")
            log.error(f"    Coloque vina.exe em: {base.resolve()}")
            log.info("\n  Gerando configs para execucao manual...")
            _generate_vina_configs(selected_nums, ligand_pdbqt, protein_dir,
                                   docking_dir, log)
            return False

    for num in selected_nums:
        target = TARGET_DATABASE[num]
        pdb_id = target["pdb"]
        tid = target["id"]

        receptor = protein_dir / f"{pdb_id}_receptor.pdbqt"
        output = docking_dir / f"result_{tid}.pdbqt"
        log_dock = docking_dir / f"log_{tid}.txt"

        if not receptor.exists():
            log.warning(f"  [AVISO] Receptor {pdb_id} nao encontrado, pulando {tid}")
            continue

        cx, cy, cz = target["grid_center"]
        sx, sy, sz = target["grid_size"]

        log.info(f"\n  [{num}] {tid}: {target['target_name']}")
        log.info(f"       Grid: center=({cx},{cy},{cz}) size=({sx},{sy},{sz})")

        try:
            if use_python_vina:
                v = Vina(sf_name='vina')
                v.set_receptor(str(receptor))
                v.set_ligand_from_file(str(ligand_pdbqt))
                v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
                v.dock(exhaustiveness=32, n_poses=10)
                energies = v.energies()
                v.write_poses(str(output), n_poses=10, overwrite=True)
                affinities = [e[0] for e in energies]
            else:
                cmd = [
                    vina_executable or "vina",
                    "--receptor", str(receptor),
                    "--ligand", str(ligand_pdbqt),
                    "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                    "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
                    "--exhaustiveness", "32",
                    "--num_modes", "10",
                    "--out", str(output),
                ]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                # Salvar stdout como log
                with open(str(log_dock), "w") as lf:
                    lf.write(result.stdout)
                affinities = _parse_vina_output(result.stdout)

            if affinities:
                best = affinities[0]
                ki = calc_ki(best)
                all_results[tid] = {
                    "number": num,
                    "pdb": pdb_id,
                    "organism": target["organism"],
                    "target_name": target["target_name"],
                    "category": target["category"],
                    "best_dG": best,
                    "ki_uM": ki,
                    "all_modes": affinities[:5],
                    "output_file": str(output),
                }
                log.info(f"       [OK] dG = {best:.1f} kcal/mol | Ki = {ki:.3f} uM")

        except Exception as e:
            log.error(f"       [ERRO] Docking falhou: {e}")

    # ----- STEP 4: Analise -----
    log.info("\n" + "=" * 60)
    log.info("STEP 4: Analise de Resultados")
    log.info("=" * 60)

    if not all_results:
        log.warning("  Nenhum resultado de docking disponivel.")
        log.info("  Verifique a preparacao dos receptores e re-execute.")
        return False

    # Tabela de resultados
    log.info(f"\n  {'Alvo':<12} {'Organismo':<18} {'Proteina':<30} "
             f"{'PDB':<6} {'dG':<10} {'Ki (uM)':<10}")
    log.info(f"  {'-'*12} {'-'*18} {'-'*30} {'-'*6} {'-'*10} {'-'*10}")

    for tid, data in sorted(all_results.items(),
                           key=lambda x: x[1]["best_dG"]):
        log.info(f"  {tid:<12} {data['organism']:<18} {data['target_name']:<30} "
                 f"{data['pdb']:<6} {data['best_dG']:<10.1f} {data['ki_uM']:<10.3f}")

    # Seletividade
    log.info(f"\n  --- Indices de Seletividade in silico ---")
    selectivity = {}
    for par_id, hum_id in SELECTIVITY_PAIRS.items():
        if par_id in all_results and hum_id in all_results:
            si = all_results[hum_id]["ki_uM"] / all_results[par_id]["ki_uM"]
            selectivity[f"{par_id}_vs_{hum_id}"] = si
            log.info(f"    SI ({par_id} vs {hum_id}): {si:.1f}x")
            if si > 1:
                log.info(f"       --> Seletivo para o parasito [OK]")
            else:
                log.info(f"       --> Sem seletividade (Ki parasito > Ki humano)")

    # Salvar tudo
    output_data = {
        "compound": COMPOUND,
        "date": datetime.now().isoformat(),
        "targets_selected": selected_nums,
        "docking_results": all_results,
        "selectivity": selectivity,
        "parameters": {
            "software": "AutoDock Vina",
            "exhaustiveness": 32,
            "num_modes": 10,
            "ligand_optimization": "MMFF94 (1000 steps CG)",
        },
    }

    json_path = results_dir / "docking_results.json"
    with open(json_path, "w") as f:
        json.dump(output_data, f, indent=2, default=str)
    log.info(f"\n  [OK] JSON: {json_path}")

    # CSV
    csv_path = results_dir / "Table_docking.csv"
    with open(csv_path, "w") as f:
        f.write("Target_ID,Organism,Protein,PDB,dG_kcal_mol,Ki_uM,Category\n")
        for tid, data in all_results.items():
            f.write(f"{tid},{data['organism']},{data['target_name']},"
                    f"{data['pdb']},{data['best_dG']:.1f},"
                    f"{data['ki_uM']:.3f},{data['category']}\n")
    log.info(f"  [OK] CSV: {csv_path}")

    # ----- STEP 5: Figuras -----
    log.info("\n" + "=" * 60)
    log.info("STEP 5: Geracao de Figuras")
    log.info("=" * 60)

    generate_figures(all_results, selectivity, figures_dir, log)
    generate_pymol_scripts(all_results, protein_dir, docking_dir, figures_dir, log)

    log.info("\n" + "=" * 60)
    log.info("[OK] PIPELINE CONCLUIDO!")
    log.info("=" * 60)
    log.info(f"\n  Resultados: {results_dir}/")
    log.info(f"  Figuras:    {figures_dir}/")
    log.info(f"  Log:        {log_file}")

    return True


def _generate_vina_configs(selected_nums, ligand, protein_dir, docking_dir, log):
    """Gera arquivos de configuracao para rodar Vina manualmente"""
    log.info("\n  Gerando configs para execucao manual do Vina:")

    for num in selected_nums:
        target = TARGET_DATABASE[num]
        tid = target["id"]
        pdb_id = target["pdb"]
        cx, cy, cz = target["grid_center"]
        sx, sy, sz = target["grid_size"]

        config_content = f"""# AutoDock Vina config - {tid} ({target['organism']})
receptor = {protein_dir}/{pdb_id}_receptor.pdbqt
ligand = {ligand}
center_x = {cx}
center_y = {cy}
center_z = {cz}
size_x = {sx}
size_y = {sy}
size_z = {sz}
exhaustiveness = 32
num_modes = 10
energy_range = 3
out = {docking_dir}/result_{tid}.pdbqt
log = {docking_dir}/log_{tid}.txt
"""
        config_path = docking_dir / f"config_{tid}.txt"
        with open(config_path, "w") as f:
            f.write(config_content)
        log.info(f"    [OK] {config_path}")
        log.info(f"      Rodar: vina --config {config_path}")


def _parse_vina_output(text):
    """Parseia stdout do Vina CLI"""
    affinities = []
    reading = False
    for line in text.splitlines():
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
    return affinities


def calc_ki(dG, R=1.987e-3, T=298.15):
    """dG (kcal/mol) --> Ki (uM)"""
    return math.exp(dG / (R * T)) * 1e6


# =============================================================================
# FIGURAS
# =============================================================================

def generate_figures(results, selectivity, fig_dir, log):
    """Gera figuras para publicacao"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import numpy as np

    # Separar parasito vs humano e ordenar
    parasite = {k: v for k, v in results.items() if v["category"] == "parasite"}
    human = {k: v for k, v in results.items() if v["category"] == "human"}

    all_sorted = list(sorted(parasite.items(), key=lambda x: x[1]["best_dG"])) + \
                 list(sorted(human.items(), key=lambda x: x[1]["best_dG"]))

    labels = []
    values = []
    colors = []

    for tid, data in all_sorted:
        labels.append(f"{data['target_name']}\n({data['organism']})")
        values.append(data["best_dG"])
        if "T. cruzi" in data["organism"]:
            colors.append("#D32F2F")
        elif "Leishmania" in data["organism"]:
            colors.append("#FF5722")
        else:
            colors.append("#1565C0")

    # ===== Figura principal: barras horizontais =====
    fig, ax = plt.subplots(figsize=(10, max(5, len(labels) * 0.8)))

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, values, color=colors, edgecolor='black',
                   linewidth=0.6, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Binding Energy (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title(f'Molecular Docking: Selectivity Profile of 2-SCN-NQ',
                 fontsize=13, fontweight='bold', pad=15)
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlim(min(values) - 1.5, 0.5)
    ax.invert_yaxis()

    # Valores nas barras
    for bar, val in zip(bars, values):
        ax.text(val - 0.15, bar.get_y() + bar.get_height() / 2.,
                f'{val:.1f}', ha='right', va='center',
                fontsize=10, fontweight='bold', color='white')

    # Legenda
    legend_elements = [
        Patch(facecolor='#D32F2F', edgecolor='black', label='T. cruzi'),
        Patch(facecolor='#FF5722', edgecolor='black', label='Leishmania spp.'),
        Patch(facecolor='#1565C0', edgecolor='black', label='H. sapiens (control)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

    # Linha separadora parasito/humano
    n_parasite = len(parasite)
    if human:
        ax.axhline(y=n_parasite - 0.5, color='gray', linestyle='--',
                   linewidth=1, alpha=0.5)
        ax.text(min(values) - 1.2, n_parasite - 0.5, 'Parasite ^ | Human v',
                fontsize=8, color='gray', va='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    for fmt in ['png', 'tiff', 'svg', 'pdf']:
        plt.savefig(fig_dir / f"Fig_selectivity.{fmt}", dpi=300,
                    bbox_inches='tight', facecolor='white')
    log.info(f"  [OK] Figura de seletividade: {fig_dir}/Fig_selectivity.*")
    plt.close()

    # ===== Figura 2: Ki comparativo (escala log) =====
    fig2, ax2 = plt.subplots(figsize=(10, max(5, len(labels) * 0.8)))

    ki_values = [results[tid]["ki_uM"] for tid, _ in all_sorted]

    bars2 = ax2.barh(y_pos, ki_values, color=colors, edgecolor='black',
                     linewidth=0.6, height=0.6)
    ax2.set_xscale('log')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_xlabel('Estimated Ki (uM) - log scale', fontsize=12, fontweight='bold')
    ax2.set_title('Estimated Inhibition Constants (Ki)',
                  fontsize=13, fontweight='bold', pad=15)
    ax2.invert_yaxis()
    ax2.legend(handles=legend_elements, loc='lower right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.tight_layout()

    for fmt in ['png', 'tiff']:
        plt.savefig(fig_dir / f"Fig_Ki_comparison.{fmt}",
                    dpi=300, bbox_inches='tight', facecolor='white')
    log.info(f"  [OK] Figura Ki: {fig_dir}/Fig_Ki_comparison.*")
    plt.close()


def generate_pymol_scripts(results, protein_dir, docking_dir, fig_dir, log):
    """Gera scripts PyMOL para cada alvo com resultado"""

    color_map = {
        "T. cruzi": "palegreen",
        "Leishmania infantum": "lightorange",
        "Leishmania major": "wheat",
        "H. sapiens": "lightblue",
    }

    for tid, data in results.items():
        pdb_id = data["pdb"]
        org = data["organism"]
        prot_color = color_map.get(org, "gray80")

        script = f"""# PyMOL - {tid} ({org})
# Uso: pymol -cq {tid}_pose.pml

load {protein_dir}/{pdb_id}_clean.pdb, receptor
load {docking_dir}/result_{tid}.pdbqt, ligand

hide everything
show cartoon, receptor
color {prot_color}, receptor
set cartoon_transparency, 0.3

show sticks, ligand
color yellow, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color orange, ligand and elem S

select site, receptor within 4 of ligand
show sticks, site
color cyan, site and elem C

show surface, site
set transparency, 0.7, site

distance hbonds, ligand, site, 3.5, mode=2
hide labels, hbonds
color yellow, hbonds

bg_color white
set ray_opaque_background, 1
set antialias, 2
orient ligand
zoom ligand, 8

ray 2400, 2400
png {fig_dir}/pose_{tid}.png, dpi=300
"""
        script_path = fig_dir / f"{tid}_pose.pml"
        with open(script_path, "w") as f:
            f.write(script)

    log.info(f"  [OK] Scripts PyMOL gerados: {fig_dir}/*_pose.pml")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Seletor interativo de alvos para docking molecular",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--targets", type=str, default=None,
                       help="Alvos separados por virgula (ex: 1,5,9)")
    parser.add_argument("--all", action="store_true",
                       help="Selecionar todos os alvos")
    parser.add_argument("--all-parasite", action="store_true",
                       help="Selecionar todos os parasitarios")
    parser.add_argument("--info", type=int, default=None,
                       help="Mostrar detalhes de um alvo (ex: --info 3)")
    parser.add_argument("--list", action="store_true",
                       help="Listar todos os alvos e sair")

    args = parser.parse_args()

    # Modo: listar alvos
    if args.list:
        display_menu()
        sys.exit(0)

    # Modo: detalhes de um alvo
    if args.info is not None:
        display_target_detail(args.info)
        sys.exit(0)

    # Modo: selecao por argumento
    if args.all:
        selected = list(TARGET_DATABASE.keys())
    elif args.all_parasite:
        selected = [n for n, d in TARGET_DATABASE.items() if d["category"] == "parasite"]
    elif args.targets:
        selected = [int(x.strip()) for x in args.targets.split(",")
                    if x.strip().isdigit() and int(x.strip()) in TARGET_DATABASE]
    else:
        # Modo interativo
        display_menu()
        selected = get_user_selection()

        if not selected:
            print("\n  Nenhum alvo selecionado. Saindo.")
            sys.exit(0)

        # Confirmar e sugerir controles
        selected = display_selection_summary(selected)

        # Confirmar execucao
        print(f"\n  Iniciar docking para {len(selected)} alvos? (S/N)")
        confirm = input("  >>> ").strip().upper()
        if confirm not in ('S', 'Y', 'SIM', 'YES', ''):
            print("\n  Cancelado.")
            sys.exit(0)

    # Executar
    print(f"\n  Alvos finais: {selected}")
    success = run_docking_for_selection(selected)

    if success:
        print("\n  [OK] Pipeline finalizado com sucesso!")
    else:
        print("\n  [AVISO] Pipeline finalizado com avisos. Verifique o log.")


if __name__ == "__main__":
    main()
