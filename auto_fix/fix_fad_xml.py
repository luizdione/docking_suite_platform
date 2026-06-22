"""
fix_fad_xml.py - Corrige cofactor_FAD_TcTR.xml deterministicamente.

Diagnostico (ver FAD_XML_SCRAMBLE_DIAGNOSIS.md):
  - 33 dos 86 atom types tem element/mass errados
  - HarmonicBondForce esta vazio (0 bonds parametrizadas)
  - Sem termo harmonico nenhuma bond e contida em dinamica -> NaN no warm-up

Estrategia de fix (sem LLM, sem internet):
  1. Para cada <Atom name="X" type="FADn">, descobrir o elemento correto
     a partir do primeiro caractere alfabetico do nome ("N9A" -> N, "C5B" -> C)
  2. Reescrever <Type element="..." mass="..."> com mass canonico daquele elemento
  3. Para cada bond da Residue, gerar <Bond class1 class2 length k> em HBF
     com length/k extraidos de uma tabela GAFF/AMBER por par de elementos
  4. Tambem gerar <Angle> e <Proper>/<Improper> minimos? Nao — bond ja resolve NaN
     na grande maioria; angle/torsion vem da template heredada se houver.
  5. Validar via simulacao no-step (criar System sem rodar steps)

Uso:
    python fix_fad_xml.py [--xml PATH] [--inplace] [--dry-run]

Default xml = docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Massa canonica e force-field defaults
ELEMENT_MASS = {
    "H": 1.00784,
    "C": 12.0108,
    "N": 14.0067,
    "O": 15.9994,
    "P": 30.9738,
    "S": 32.065,
}

# Parametros de bond por par de elementos (GAFF2 typical, OpenMM units: nm + kJ/mol/nm^2)
# Fontes: GAFF2 paper, OpenFF 2.x SMIRNOFF defaults. Valores conservadores para
# estabilizar geometria — NAO sao parametros "publicaveis" (faltam refit), mas
# bastam para nao explodir NaN no warm-up.
BOND_PARAMS = {  # (elem1, elem2) -> (length_nm, k_kJ/mol/nm^2)
    ("C", "C"): (0.1525, 259408.0),   # single C-C
    ("C", "N"): (0.1470, 282001.6),   # single C-N
    ("C", "O"): (0.1410, 267776.0),   # single C-O
    ("C", "H"): (0.1090, 284512.0),
    ("N", "H"): (0.1010, 363171.2),
    ("O", "H"): (0.0960, 462750.4),
    ("C", "P"): (0.1830, 188280.0),
    ("O", "P"): (0.1610, 188280.0),
    ("N", "N"): (0.1410, 251040.0),
    ("N", "O"): (0.1410, 251040.0),
    ("P", "O"): (0.1610, 188280.0),
    ("P", "H"): (0.1420, 188280.0),
    ("S", "C"): (0.1810, 167360.0),
    ("S", "H"): (0.1336, 230120.0),
    ("S", "S"): (0.2038, 138072.0),
    ("S", "O"): (0.1810, 188280.0),
    # Aromatic / ring overrides (when in adenine/isoalloxazine rings):
    # detected by name pattern N1A/C2A/N3A/C4A/C5A/C6A/C8A/N7A/N9A (adenine)
    # or N1/C2/N3/C4/C4X/C5X/N5/C6/C7/C8/C9/C9A/N10/C10A (isoalloxazine)
    # We don't differentiate here — bonds inside these rings get the same
    # values; for higher accuracy use openff-toolkit or antechamber.
}

DEFAULT_XML = Path(__file__).resolve().parents[1] / "docking_project" / \
    "figures" / "smoke_test_fad" / "cofactor_FAD_TcTR.xml"


def correct_element(atom_name: str) -> str:
    """First alphabetic char of the atom name, upper-cased."""
    for ch in atom_name:
        if ch.isalpha():
            return ch.upper()
    raise ValueError(f"No alpha char in atom name {atom_name!r}")


def parse_xml(xml_text: str) -> dict:
    """Lightweight regex-based parser — XML lib breaks on the 1-line collapsed format."""
    info: dict = {}

    # Atom name -> type
    info["atom_to_type"] = {
        m.group(1): m.group(2)
        for m in re.finditer(r'<Atom name="([^"]+)" type="(FAD\d+)"', xml_text)
    }

    # Type defs: name, element, mass
    types: Dict[str, Tuple[str, float, str]] = {}
    for m in re.finditer(
        r'<Type name="(FAD\d+)" element="([^"]+)" mass="([^"]+)" class="(FAD\d+)"\s*/>',
        xml_text,
    ):
        types[m.group(1)] = (m.group(2), float(m.group(3)), m.group(4))
    info["types"] = types

    # Residue bonds (atom_name pairs)
    res_match = re.search(r'<Residue name="FAD">(.*?)</Residue>', xml_text, re.DOTALL)
    if res_match:
        info["res_block"] = res_match.group(0)
        info["res_bonds"] = [
            (m.group(1), m.group(2))
            for m in re.finditer(
                r'<Bond atomName1="([^"]+)" atomName2="([^"]+)"',
                res_match.group(1),
            )
        ]
    else:
        info["res_block"] = ""
        info["res_bonds"] = []

    # Existing HBF entries
    hbf = re.search(r'<HarmonicBondForce>(.*?)</HarmonicBondForce>', xml_text, re.DOTALL)
    info["hbf_block"] = hbf.group(0) if hbf else ""
    info["hbf_bonds"] = []
    if hbf:
        # FIX 2026-05-21: o XML do FAD grava bonds como type1/type2 (nao
        # class1/class2). O regex antigo so casava class1/class2 -> lia 0 das 91
        # bonds existentes e o builder injetava 91 DUPLICATAS (dobrando a mola).
        # Agora aceita os dois formatos.
        for m in re.finditer(
            r'<Bond (?:class1|type1)="(FAD\d+)" (?:class2|type2)="(FAD\d+)" '
            r'length="([^"]+)" k="([^"]+)"',
            hbf.group(1),
        ):
            info["hbf_bonds"].append((m.group(1), m.group(2),
                                      float(m.group(3)), float(m.group(4))))
    return info


def patch_atom_types(xml_text: str, atom_to_type, type_to_corrected) -> Tuple[str, int]:
    """Rewrite <Type ... element="X" mass="Y" ...> based on corrected element."""
    n_changed = 0

    def repl(m: re.Match) -> str:
        nonlocal n_changed
        tp = m.group(1)
        old_elem = m.group(2)
        old_mass = float(m.group(3))
        klass = m.group(4)
        corr = type_to_corrected.get(tp)
        if corr is None:
            return m.group(0)  # unchanged (type not used by any atom)
        new_elem, new_mass = corr
        if new_elem != old_elem or abs(new_mass - old_mass) > 0.001:
            n_changed += 1
            return (f'<Type name="{tp}" element="{new_elem}" '
                    f'mass="{new_mass:.4f}" class="{klass}" />')
        return m.group(0)

    new_text = re.sub(
        r'<Type name="(FAD\d+)" element="([^"]+)" mass="([^"]+)" class="(FAD\d+)"\s*/>',
        repl,
        xml_text,
    )
    return new_text, n_changed


def bond_params_for(e1: str, e2: str) -> Tuple[float, float] | None:
    """Lookup (length, k) for a pair of elements. Returns None if unsupported."""
    key = tuple(sorted([e1, e2]))
    rev = (key[1], key[0])
    return BOND_PARAMS.get(key) or BOND_PARAMS.get(rev)


def build_missing_hbf_entries(info: dict, atom_correct_element) -> List[str]:
    """For each residue bond not present in HBF, emit a <Bond ...> XML entry."""
    # SEGURANCA 2026-05-21: se o HBF ja tem bonds (o FAD tem 91 validas), NAO
    # injetar nada. As bonds NAO estao faltando -- o diagnostico de 14/05 estava
    # errado. Injetar duplicaria a mola harmonica. Ver BRIEFING_2026-05-21.
    if info["hbf_bonds"]:
        print(f"  [skip] HBF ja tem {len(info['hbf_bonds'])} bonds; "
              f"nenhuma injetada (evita duplicatas).")
        return []
    have = {frozenset({a, b}) for a, b, _, _ in info["hbf_bonds"]}
    new_lines: List[str] = []
    for a, b in info["res_bonds"]:
        ta = info["atom_to_type"].get(a)
        tb = info["atom_to_type"].get(b)
        if not (ta and tb):
            continue
        if frozenset({ta, tb}) in have:
            continue
        ea = atom_correct_element(a)
        eb = atom_correct_element(b)
        params = bond_params_for(ea, eb)
        if params is None:
            # Unknown pair — fall back to generic stiff bond at 0.15 nm
            length, k = 0.15, 250000.0
        else:
            length, k = params
        new_lines.append(
            f'  <Bond class1="{ta}" class2="{tb}" length="{length:.4f}" '
            f'k="{k:.1f}"/>'
        )
        have.add(frozenset({ta, tb}))
    return new_lines


def patch_hbf(xml_text: str, new_entries: List[str]) -> Tuple[str, int]:
    """Inject new <Bond> entries into <HarmonicBondForce>."""
    if not new_entries:
        return xml_text, 0
    injected = "\n" + "\n".join(new_entries) + "\n"
    # Handle both <HarmonicBondForce/> and <HarmonicBondForce>...</HarmonicBondForce>
    if "<HarmonicBondForce/>" in xml_text:
        new_text = xml_text.replace(
            "<HarmonicBondForce/>",
            "<HarmonicBondForce>" + injected + "</HarmonicBondForce>",
        )
        return new_text, len(new_entries)
    if "<HarmonicBondForce></HarmonicBondForce>" in xml_text:
        new_text = xml_text.replace(
            "<HarmonicBondForce></HarmonicBondForce>",
            "<HarmonicBondForce>" + injected + "</HarmonicBondForce>",
        )
        return new_text, len(new_entries)
    # Else: insert before the closing tag
    if "</HarmonicBondForce>" in xml_text:
        new_text = xml_text.replace(
            "</HarmonicBondForce>",
            injected + "</HarmonicBondForce>",
        )
        return new_text, len(new_entries)
    # No HBF block at all — create one before </ForceField>
    block = ("<HarmonicBondForce>" + injected + "</HarmonicBondForce>")
    new_text = xml_text.replace("</ForceField>", block + "</ForceField>")
    return new_text, len(new_entries)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xml", default=str(DEFAULT_XML),
                    help=f"Path to FAD XML (default: {DEFAULT_XML.name})")
    ap.add_argument("--inplace", action="store_true",
                    help="Overwrite the input XML (backup with .bak suffix)")
    ap.add_argument("--out", default=None,
                    help="Output XML (default: same dir + _fixed.xml)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Just report what would change, do not write")
    args = ap.parse_args(argv)

    print("=" * 72)
    print("AVISO (2026-05-21): este script e' um STOPGAP e NAO e' o fix definitivo.")
    print("  - Corrige elemento/massa dos 33 atom types errados.")
    print("  - NAO corrige LJ (sigma/epsilon), cargas nem torsions (todos ainda")
    print("    parametrizados para o atomo ERRADO / vazios).")
    print("  - Corrigir SO o XML desincroniza da coluna element do PDB e quebra o")
    print("    template matching ('No template found for residue FAD').")
    print("  => Fix correto: rode  python regen_fad_xml.py  (rota SMIRNOFF completa).")
    print("=" * 72)

    xml_path = Path(args.xml)
    if not xml_path.exists():
        print(f"FAIL: XML not found: {xml_path}", file=sys.stderr)
        return 2

    xml_text = xml_path.read_text(encoding="utf-8")
    info = parse_xml(xml_text)

    print(f"XML: {xml_path}")
    print(f"  Atoms: {len(info['atom_to_type'])}")
    print(f"  Types: {len(info['types'])}")
    print(f"  Residue bonds: {len(info['res_bonds'])}")
    print(f"  HBF bonds (existing): {len(info['hbf_bonds'])}")

    # Build atom-name -> correct-element map
    atom_to_elem = {a: correct_element(a) for a in info["atom_to_type"]}
    # Build type -> (correct_element, correct_mass) using atom_to_type
    type_to_corrected: Dict[str, Tuple[str, float]] = {}
    n_already_ok = 0
    for atom, tp in info["atom_to_type"].items():
        elem = atom_to_elem[atom]
        mass = ELEMENT_MASS.get(elem, 0.0)
        type_to_corrected[tp] = (elem, mass)
        old = info["types"].get(tp)
        if old and old[0] == elem and abs(old[1] - mass) < 0.001:
            n_already_ok += 1
    print(f"  Types already correct: {n_already_ok}/{len(info['types'])}")
    print(f"  Types to fix: {len(info['types']) - n_already_ok}")

    new_text, n_changed = patch_atom_types(xml_text, info["atom_to_type"],
                                            type_to_corrected)
    print(f"  Type entries rewritten: {n_changed}")

    # After fixing types, the type->element map is stable; build HBF entries
    info_after = parse_xml(new_text)
    new_entries = build_missing_hbf_entries(info_after, lambda a: atom_to_elem[a])
    print(f"  Missing HBF bonds to add: {len(new_entries)}")

    new_text, n_hbf = patch_hbf(new_text, new_entries)
    print(f"  HBF entries injected: {n_hbf}")

    if args.dry_run:
        print("\n[--dry-run] No file written.")
        return 0

    out = (xml_path if args.inplace
           else (Path(args.out) if args.out else
                 xml_path.with_name(xml_path.stem + "_fixed.xml")))
    if args.inplace:
        backup = xml_path.with_suffix(xml_path.suffix + ".bak")
        if not backup.exists():
            shutil.copy2(xml_path, backup)
            print(f"  Backup saved: {backup.name}")
    out.write_text(new_text, encoding="utf-8")
    print(f"  Written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
