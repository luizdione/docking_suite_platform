#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
regen_ligand_xml.py  (v1-2026-05-25)

Regenera o force-field E o template de QUALQUER LIGANTE para a MD, da mesma forma
robusta que o regen_fad_xml.py faz para o cofator: separa POSE (vem do docking) de
QUIMICA (ordens de ligacao vem de uma fonte real). Resolve o bug que estava no
lig_TcTR_smirnoff.xml: 0 torsions + ligante super-protonado.

CAUSA-RAIZ (descoberta em 2026-05-25): o ligand_TcTR.sdf / lig_TcTR.pdb usados na
MD tinham TODAS as 37 ligacoes como ordem 1 (single) -- ordens de ligacao perdidas.
Sem duplas/aromaticos, o RDKit AddHs super-protonou o esqueleto C11NO2S para 21 H
(C11H21NO2S, 36 atomos) -- uma "naftoquinona" saturada quimicamente impossivel. O
SMIRNOFF entao gera um FF sem rigidez de anel (poucos/zero torsions proprios) -> a
MD roda sobre a molecula ERRADA.

A POSE (esqueleto pesado: 15 atomos C11NO2S no referencial da proteina) esta
CORRETA no lig_TcTR.pdb -- so as ordens e os H estavam errados. Logo:
  1. Le os atomos PESADOS da pose docada (lig_TcTR.pdb), descartando os H errados.
  2. Percebe a conectividade pesada por DISTANCIA (geometria real -> grafo certo).
  3. Atribui ORDENS DE LIGACAO a partir de uma fonte real (prioridade):
       --smiles  >  --sdf (se tiver ordens reais)  >  --mol2  >  --ccd  >  default
     via AssignBondOrdersFromTemplate (VF2). Para o 2-SCN-NQ o default ja e o SMILES
     verificado a partir do .mol2 do docking.
  4. Re-adiciona os H corretos (AddHs com coordenadas) -> molecula certa na pose.
  5. Parametriza via SMIRNOFF (reusa build_xml do regen_fad_xml -> elementos/massas/
     LJ/cargas/bonds/angles/TORSIONS), com prefixo de tipo = nome do residuo (LIG)
     para NAO colidir com os tipos do FAD quando os dois XML carregam juntos.
  6. Grava lig_<tag>_smirnoff.xml + lig_<tag>.pdb (com CONECT correto) e auto-valida
     com a forcefield COMPLETA (amber14 + tip3p + FAD_v3 + este ligante).

Para um COMPOSTO NOVO contra a TcTR: rode o docking, salve o .mol2 (ou tenha o
SMILES) e rode:
    python regen_ligand_xml.py --pdb <pose_docada.pdb> --mol2 <novo.mol2> \
        --name LIG --tag <novo_tag>
Depois passe --lig-xml lig_<novo_tag>_smirnoff.xml para o md_clean.py.

Rode no env `datascience`:  python regen_ligand_xml.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Reusa a rota SMIRNOFF->XML JA VALIDADA do regen_fad_xml (fonte unica de verdade).
# build_xml usa o nome do residuo como prefixo de tipo (FAD0.. / LIG0..) -> sem colisao.
try:
    from regen_fad_xml import build_xml, _parse_ccd_smiles
except Exception as e:  # pragma: no cover
    sys.exit(f"FALHA ao importar de regen_fad_xml.py (deve estar na mesma pasta): {e}")

SCRIPT_VERSION = "v1-2026-05-25"
HERE = Path(__file__).resolve().parent
SMOKE = HERE / "docking_project" / "figures" / "smoke_test_fad"

DEF_PDB = SMOKE / "lig_TcTR.pdb"                                   # pose docada (heavy ok)
DEF_MOL2 = HERE / "docking_project" / "ligand" / "2-SCN-NQ" / "2-SCN-NQ.mol2"
FAD_V3_XML = SMOKE / "cofactor_FAD_TcTR_v3.xml"

HH_CUT = 1.95   # A, corte ligacao pesada-pesada
# SMILES verificado do 2-SCN-NQ (2-tiocianato-1,4-naftoquinona), derivado do .mol2
# do docking (C=O x2, anel aromatico, S-C#N nitrila). C11H5NO2S.
DEFAULT_SMILES = "O=C1c2ccccc2C(=O)C(SC#N)=C1"


# ─────────────────────────────────────────────────────────────────────────────
# Fontes de ORDENS DE LIGACAO (SMILES de referencia)
# ─────────────────────────────────────────────────────────────────────────────
def smiles_from_mol2(path: Path):
    """Deriva SMILES de um .mol2 (artefato natural do docking) via RDKit."""
    from rdkit import Chem
    for san in (True, False):
        try:
            m = Chem.MolFromMol2File(str(path), sanitize=san, removeHs=True)
            if m is not None:
                if not san:
                    Chem.SanitizeMol(m)
                smi = Chem.MolToSmiles(Chem.RemoveHs(m))
                print(f"  [mol2] SMILES derivado de {path.name}: {smi}")
                return smi
        except Exception as e:
            last = e
    print(f"  [mol2] nao consegui ler {path.name} ({str(last)[:90] if 'last' in dir() else ''})")
    return None


def smiles_from_sdf(path: Path):
    """Deriva SMILES de um SDF, MAS recusa se todas as ligacoes forem single
    (sinal de ordens perdidas -- o bug que estamos corrigindo)."""
    from rdkit import Chem
    # checagem barata antes de confiar: bloco de contagem + ordens
    try:
        lines = path.read_text(errors="replace").splitlines()
        counts = lines[3].split()
        nbond = int(counts[1])
        orders = []
        for ln in lines[4 + int(counts[0]): 4 + int(counts[0]) + nbond]:
            orders.append(ln[6:9].strip())
        if orders and all(o == "1" for o in orders):
            print(f"  [sdf] {path.name}: TODAS as {nbond} ligacoes sao single "
                  f"(ordens perdidas) -> RECUSADO. Use --mol2 ou --smiles.")
            return None
    except Exception:
        pass
    m = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=True)
    if m is None:
        print(f"  [sdf] RDKit nao leu {path.name}")
        return None
    smi = Chem.MolToSmiles(m)
    print(f"  [sdf] SMILES de {path.name}: {smi}")
    return smi


def fetch_ccd_smiles(code: str):
    """Baixa o SMILES canonico do CCD do RCSB para um codigo de 3 letras."""
    import urllib.request
    for tmpl in ("https://files.rcsb.org/ligands/download/{}.cif",
                 "https://files.rcsb.org/ligands/view/{}.cif"):
        url = tmpl.format(code.upper())
        try:
            print(f"  [CCD] baixando {url}")
            with urllib.request.urlopen(url, timeout=30) as r:
                smi = _parse_ccd_smiles(r.read().decode("utf-8", "replace"))
                if smi:
                    return smi
        except Exception as e:
            last = e
    print(f"  [CCD] falhou para {code}: {str(last)[:90] if 'last' in dir() else ''}")
    return None


def resolve_smiles(args):
    if args.smiles:
        print(f"  [fonte] SMILES explicito: {args.smiles}")
        return args.smiles
    if args.sdf and Path(args.sdf).exists():
        s = smiles_from_sdf(Path(args.sdf))
        if s:
            return s
    if args.mol2 and Path(args.mol2).exists():
        s = smiles_from_mol2(Path(args.mol2))
        if s:
            return s
    if args.ccd:
        s = fetch_ccd_smiles(args.ccd)
        if s:
            return s
    if args.allow_default:
        print(f"  [fonte] usando SMILES DEFAULT do 2-SCN-NQ: {DEFAULT_SMILES}")
        print("          (para outro composto passe --smiles/--mol2/--sdf/--ccd)")
        return DEFAULT_SMILES
    sys.exit("Sem fonte de ordens de ligacao. Passe --smiles, --mol2, --sdf ou --ccd "
             "(ou --allow-default para o 2-SCN-NQ).")


# ─────────────────────────────────────────────────────────────────────────────
# Construir o ligante: pose pesada (docking) + ordens (SMILES) + H corretos
# ─────────────────────────────────────────────────────────────────────────────
def read_heavy_pdb(pdb: Path):
    """Le SO os atomos pesados (descarta H, mesmo os super-protonados errados)."""
    names, elems, xyz = [], [], []
    for ln in pdb.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        nm = ln[12:16].strip()
        el = (ln[76:78].strip() or "".join(c for c in nm if c.isalpha())[:1]).capitalize()
        if el == "H":
            continue
        names.append(nm)
        elems.append(el)
        xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return names, elems, xyz


def build_ligand(pdb: Path, smiles: str):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    from rdkit.Geometry import Point3D

    xn, xe, xyz = read_heavy_pdb(pdb)
    xyz = np.array(xyz)
    nH = len(xn)
    print(f"  [pose] {pdb.name}: {nH} atomos pesados (H descartados)")
    if nH < 3:
        raise RuntimeError(f"poucos atomos pesados ({nH}) em {pdb.name}")

    # conectividade pesada por distancia (geometria real da pose docada)
    D = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
    hbonds = [(i, j) for i in range(nH) for j in range(i + 1, nH) if 0.4 < D[i, j] < HH_CUT]
    print(f"  [pose] {len(hbonds)} ligacoes pesadas percebidas "
          f"({len(hbonds) - nH + 1} aneis)")

    rw = Chem.RWMol()
    for e in xe:
        rw.AddAtom(Chem.Atom(e))
    for i, j in hbonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(nH)
    for i, (x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)

    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        raise RuntimeError(f"SMILES invalido: {smiles}")
    nref_heavy = ref.GetNumHeavyAtoms()
    if nref_heavy != nH:
        raise RuntimeError(
            f"SMILES tem {nref_heavy} atomos pesados, mas a pose tem {nH}. "
            "Confirme que o SMILES e do MESMO composto da pose docada.")
    mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)  # VF2: casa o grafo pesado
    mol.UpdatePropertyCache(strict=False)

    molH = Chem.AddHs(mol, addCoords=True)
    Chem.SanitizeMol(molH)

    names = list(xn)
    hc = 0
    for a in molH.GetAtoms():
        if a.GetIdx() >= nH:
            hc += 1
            names.append(f"H{hc}")
    formula = rdMolDescriptors.CalcMolFormula(molH)
    print(f"  [mol] formula={formula} | {molH.GetNumAtoms()} atomos "
          f"({nH} pesados + {hc} H) | carga formal={Chem.GetFormalCharge(molH)}")
    return molH, names


def write_pdb(molH, names, residue: str, out_pdb: Path):
    conf = molH.GetConformer()
    lines = []
    for i, a in enumerate(molH.GetAtoms()):
        p = conf.GetAtomPosition(i)
        el = a.GetSymbol()
        nm = names[i]
        namecol = f" {nm:<3}" if len(nm) < 4 else nm[:4]
        lines.append(f"HETATM{i+1:5d} {namecol} {residue:<3} A   1    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2}")
    nbrs = {i: [] for i in range(molH.GetNumAtoms())}
    for b in molH.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        nbrs[a1].append(a2); nbrs[a2].append(a1)
    for i in range(molH.GetNumAtoms()):
        if nbrs[i]:
            lines.append("CONECT" + f"{i+1:5d}" + "".join(f"{j+1:5d}" for j in sorted(nbrs[i])))
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n")
    print(f"  [pdb] {out_pdb.name} gravado ({molH.GetNumAtoms()} atomos, "
          f"{sum(len(v) for v in nbrs.values())//2} ligacoes em CONECT)")


# ─────────────────────────────────────────────────────────────────────────────
def self_validate(out_xml: Path, out_pdb: Path, residue: str):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    print("\n=== AUTO-VALIDACAO (forcefield COMPLETA, igual a MD) ===")
    try:
        ff_files = ["amber14-all.xml", "amber14/tip3p.xml"]
        if FAD_V3_XML.exists():
            ff_files.append(str(FAD_V3_XML))
        ff_files.append(str(out_xml))
        ff = app.ForceField(*ff_files)
        pdb = app.PDBFile(str(out_pdb))
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=None, rigidWater=False)
    except Exception as e:
        print(f"  FALHA createSystem (mesmo erro que a MD veria): {str(e)[:200]}")
        return False
    ntor = sum(f.getNumTorsions() for f in system.getForces()
               if isinstance(f, mm.PeriodicTorsionForce))
    nbond = sum(f.getNumBonds() for f in system.getForces()
                if isinstance(f, mm.HarmonicBondForce))
    print(f"  createSystem OK | bonds={nbond} | torsions={ntor}")
    integ = mm.LangevinMiddleIntegrator(5 * unit.kelvin, 1 / unit.picosecond,
                                        0.00025 * unit.picoseconds)
    sim = app.Simulation(pdb.topology, system, integ, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=300)
    e = np.nan
    for _ in range(8):
        sim.step(50)
        e = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        if not np.isfinite(e):
            print("  WARM-UP ligante isolado: NaN"); return False
    ok = ntor > 0
    print(f"  WARM-UP ligante (400 passos @5K): OK (E={e:.3e} kJ/mol)")
    print("  RESULTADO:", "PASS" if ok else "ATENCAO (0 torsions -> revise)")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", type=Path, default=DEF_PDB,
                    help="pose docada do ligante (fonte da GEOMETRIA pesada)")
    ap.add_argument("--smiles", default=None, help="SMILES com ordens corretas (preferido p/ composto novo)")
    ap.add_argument("--mol2", type=Path, default=DEF_MOL2, help="MOL2 do docking (deriva SMILES)")
    ap.add_argument("--sdf", type=Path, default=None, help="SDF com ordens reais (deriva SMILES)")
    ap.add_argument("--ccd", default=None, help="codigo CCD de 3 letras (RCSB)")
    ap.add_argument("--name", default="LIG", help="nome do residuo (default LIG)")
    ap.add_argument("--tag", default="TcTR", help="sufixo de saida: lig_<tag>_smirnoff.xml")
    ap.add_argument("--out-dir", type=Path, default=SMOKE)
    ap.add_argument("--allow-default", action="store_true",
                    help="se nenhuma fonte funcionar, usar o SMILES default do 2-SCN-NQ")
    ap.add_argument("--no-sync", action="store_true",
                    help="NAO sobrescrever lig_<tag>_smirnoff.xml original (so grava *_v2)")
    a = ap.parse_args()

    # por conveniencia do projeto TcTR, habilita o default automaticamente quando
    # nenhuma fonte explicita foi passada e o pdb e o do 2-SCN-NQ
    if not (a.smiles or a.sdf or a.ccd) and not a.mol2.exists():
        a.allow_default = True

    print(f"[regen_ligand_xml.py {SCRIPT_VERSION}]  (pose docking + ordens reais -> SMIRNOFF)")
    if not a.pdb.exists():
        sys.exit(f"FALTA a pose: {a.pdb}")
    out_xml = a.out_dir / f"lig_{a.tag}_smirnoff_v2.xml"
    out_pdb = a.out_dir / f"lig_{a.tag}_v2.pdb"
    print(f"Pose  : {a.pdb}")
    print(f"Saidas: {out_xml.name} + {out_pdb.name}\n")

    smiles = resolve_smiles(a)
    molH, names = build_ligand(a.pdb, smiles)
    write_pdb(molH, names, a.name, out_pdb)
    xml = build_xml(molH, names, residue_name=a.name)
    out_xml.write_text(xml, encoding="utf-8")
    print(f"\nGravado: {out_xml.name}")

    ok = self_validate(out_xml, out_pdb, a.name)

    # Sincroniza com os arquivos que o md_clean usa por default (faz backup 1x),
    # a menos que --no-sync. Assim o md_clean passa a usar o ligante CORRETO.
    if ok and not a.no_sync:
        import shutil
        for _new, _orig in [(out_xml, a.out_dir / f"lig_{a.tag}_smirnoff.xml"),
                            (out_pdb, a.out_dir / f"lig_{a.tag}.pdb")]:
            _bak = _orig.with_suffix(_orig.suffix + ".broken_bak")
            if _orig.exists() and not _bak.exists():
                shutil.copy2(_orig, _bak)
            shutil.copy2(_new, _orig)
            print(f"  [sync] {_orig.name} <- {_new.name}")
    elif not ok:
        print("  [sync] PULADO (validacao falhou; nao sobrescrevo os originais).")

    print("\n" + "=" * 64)
    print("PROXIMO PASSO:")
    print(f"  python validate_ff.py {out_xml} --pdb {out_pdb}")
    print("  python md_simulation.py --fresh --ns 0.2 --tag prod_v4 --platform CUDA \\")
    print(f"      --fad-xml {FAD_V3_XML} --lig-xml lig_{a.tag}_smirnoff.xml")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
