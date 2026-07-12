#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
regen_fad_xml.py  (vXTAL-2026-05-22)

Regenera o force-field E o template do FAD a partir da GEOMETRIA
CRISTALOGRAFICA (cofactor_FAD_0_xtal.pdb) -- a unica fonte CORRETA.

Causa-raiz (descoberta em 2026-05-22): o template antigo
`cofactor_FAD_TcTR.pdb` (e o XML derivado dele) tem coordenadas, nomes e CONECT
EMBARALHADOS -- fosforo com 1 ligacao, ligacoes O-O, metila com 3 vizinhos, etc.
(Kabsch RMSD template->xtal = 8.7 A). Por isso TODA tentativa de regenerar a
partir do template falhou ("No matching found"). O cristal, ao contrario, esta
intacto: PA->4 O, P->4 O, 6 aneis, 0 erros quimicos.

Estrategia desta versao:
  1. Le os 53 atomos pesados do XTAL (nomes + elementos + coordenadas corretos).
  2. Percebe a conectividade pesada por DISTANCIA (geometria real -> grafo certo).
  3. Atribui ordens de ligacao com o SMILES canonico do CCD do PDB
     (AssignBondOrdersFromTemplate -- agora casa, pois a conectividade esta certa).
  4. Adiciona os H (RDKit AddHs, com coordenadas) -> 86 atomos, FAD neutro.
  5. Parametriza via SMIRNOFF (mesma rota do ligante) -> elementos/massas/LJ/
     cargas/bonds/angles/TORSIONS corretos.
  6. Grava:
       cofactor_FAD_TcTR_v3.pdb  (86 atomos, CONECT CORRETO -> conserta tambem a
                                  colocacao de H no md_clean.build_canonical_fad_topology)
       cofactor_FAD_TcTR_v3.xml  (ForceField OpenMM, drop-in)
  7. Auto-valida (createSystem + warm-up curto do FAD isolado).

Rode no env `datascience`:  python regen_fad_xml.py
"""
from __future__ import annotations

import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPT_VERSION = "vNAGL-2026-05-23"  # se nao aparecer no log, OneDrive nao sincronizou

HERE = Path(__file__).resolve().parent
SMOKE = HERE / "docking_project" / "figures" / "smoke_test_fad"
XTAL_PDB = SMOKE / "cofactor_FAD_0_xtal.pdb"      # 53 heavy CRISTALOGRAFICOS (fonte correta)
OUT_PDB = SMOKE / "cofactor_FAD_TcTR_v3.pdb"      # template 86 atomos (regerado)
OUT_XML = SMOKE / "cofactor_FAD_TcTR_v3.xml"      # force-field drop-in

CCD_URLS = ["https://files.rcsb.org/ligands/download/FAD.cif",
            "https://files.rcsb.org/ligands/view/FAD.cif"]
LOCAL_CCD = [HERE / "FAD.cif", SMOKE / "FAD.cif"]

EXPECTED = {"C": 27, "H": 33, "N": 9, "O": 15, "P": 2}  # FAD neutro
KCAL_KJ = 4.184
ANG_NM = 0.1
HH_CUT = 1.95   # A, corte ligacao pesada-pesada
XH_CUT = 1.30   # A, corte X-H


# ─────────────────────────────────────────────────────────────────────────────
# CCD: SMILES canonico do FAD (para ordens de ligacao)
# ─────────────────────────────────────────────────────────────────────────────
def _parse_ccd_smiles(text: str):
    lines = text.splitlines()

    def read_loop(prefix):
        for i, ln in enumerate(lines):
            if ln.strip() == "loop_":
                j = i + 1
                headers = []
                while j < len(lines) and lines[j].strip().startswith("_"):
                    headers.append(lines[j].strip()); j += 1
                if headers and headers[0].startswith(prefix):
                    cols = {h.split(".", 1)[1]: k for k, h in enumerate(headers)}
                    rows = []
                    while j < len(lines):
                        s = lines[j].strip()
                        if s in ("", "#") or s.startswith(("loop_", "_")):
                            break
                        toks, cur, q = [], "", None
                        for ch in s:
                            if q:
                                if ch == q: q = None
                                else: cur += ch
                            elif ch in "'\"": q = ch
                            elif ch == " ":
                                if cur: toks.append(cur); cur = ""
                            else: cur += ch
                        if cur: toks.append(cur)
                        if toks: rows.append(toks)
                        j += 1
                    return cols, rows
        return None, None

    cols, rows = read_loop("_pdbx_chem_comp_descriptor.")
    if cols and rows and "type" in cols and "descriptor" in cols:
        ti, di = cols["type"], cols["descriptor"]
        smis = [(r[ti].upper(), r[di]) for r in rows if len(r) > max(ti, di)]
        for want in ("SMILES_CANONICAL", "SMILES"):
            for t, d in smis:
                if want in t and d:
                    return d
    return None


def load_ccd_smiles():
    import urllib.request
    for p in LOCAL_CCD:
        if p.exists():
            print(f"  [CCD] arquivo local {p.name}")
            return _parse_ccd_smiles(p.read_text())
    last = None
    for url in CCD_URLS:
        try:
            print(f"  [CCD] baixando {url}")
            with urllib.request.urlopen(url, timeout=30) as r:
                return _parse_ccd_smiles(r.read().decode("utf-8", "replace"))
        except Exception as e:
            last = e
    raise RuntimeError(f"Nao obtive o CCD do FAD (sem internet?). Baixe FAD.cif "
                       f"e ponha ao lado do script. Ultimo erro: {last}")


# ─────────────────────────────────────────────────────────────────────────────
# Construir o FAD a partir do XTAL (geometria correta) + ordens do CCD + H
# ─────────────────────────────────────────────────────────────────────────────
def build_fad_from_xtal(smiles):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    # 1) ler 53 atomos pesados do xtal
    xn, xe, xyz = [], [], []
    for ln in XTAL_PDB.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            nm = ln[12:16].strip()
            el = (ln[76:78].strip() or "".join(c for c in nm if c.isalpha())[:1]).capitalize()
            xn.append(nm); xe.append(el)
            xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    xyz = np.array(xyz); nH = len(xn)
    if nH != 53:
        raise RuntimeError(f"esperava 53 heavy no xtal, achei {nH} em {XTAL_PDB.name}")

    # 2) conectividade pesada por distancia (geometria real)
    D = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
    hbonds = [(i, j) for i in range(nH) for j in range(i + 1, nH) if 0.4 < D[i, j] < HH_CUT]
    print(f"  [xtal] {nH} heavy, {len(hbonds)} ligacoes pesadas percebidas, "
          f"{len(hbonds) - nH + 1} aneis")

    # 3) RWMol pesado + conformacao
    rw = Chem.RWMol()
    for e in xe:
        rw.AddAtom(Chem.Atom(e))
    for i, j in hbonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(nH)
    for i, (x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf, assignId=True)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)

    # 4) ordens de ligacao via SMILES do CCD (agora casa: conectividade correta)
    ref = Chem.MolFromSmiles(smiles) if smiles else None
    if ref is None:
        raise RuntimeError(f"SMILES do CCD invalido: {smiles}")
    mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)

    # 5) FAD neutro + adicionar H com coordenadas
    for a in mol.GetAtoms():
        a.SetFormalCharge(0)
    mol.UpdatePropertyCache(strict=False)
    molH = Chem.AddHs(mol, addCoords=True)
    Chem.SanitizeMol(molH)

    # 6) nomes: 53 pesados = nomes do xtal ; H = H1..Hk
    names = list(xn)
    hc = 0
    for a in molH.GetAtoms():
        if a.GetIdx() >= nH:
            hc += 1
            names.append(f"H{hc}")
    from rdkit.Chem import rdMolDescriptors
    formula = rdMolDescriptors.CalcMolFormula(molH)
    print(f"  [mol] formula={formula} (esperado C27H33N9O15P2) | {molH.GetNumAtoms()} atomos")
    if molH.GetNumAtoms() != 86:
        print(f"  [mol] AVISO: {molH.GetNumAtoms()} atomos (esperado 86). Revise protonacao.")
    return molH, names


# ─────────────────────────────────────────────────────────────────────────────
# SMIRNOFF labels -> OpenMM ForceField XML  (mesma rota do _build_lig_xml_smirnoff)
# ─────────────────────────────────────────────────────────────────────────────
def build_xml(mol_rdkit, names, residue_name="FAD"):
    import numpy as np
    try:
        from openff.toolkit import (Molecule as OFFMol, ForceField as OFFff,
                                    Topology as OFFTop)
    except OSError as e:
        raise RuntimeError(
            "Falha ao importar openff.toolkit -- quase certamente o PyTorch do "
            "openff-nagl esta com DLL quebrada (fbgemm.dll, WinError 127). "
            "CONSERTE: instale o Visual C++ Redistributable x64 (https://aka.ms/vs/17/"
            "release/vc_redist.x64.exe) e teste 'python -c \"import torch\"'. "
            "OU desinstale o NAGL p/ destravar com Gasteiger: conda remove openff-nagl "
            "openff-nagl-models pytorch pytorch-lightning torchmetrics libtorch. "
            f"[orig: {str(e)[:90]}]") from e
    from openff.units import unit as ou

    mol_off = OFFMol.from_rdkit(mol_rdkit, allow_undefined_stereo=True)
    n = mol_off.n_atoms
    if n != len(names):
        raise RuntimeError(f"n_atoms ({n}) != nomes ({len(names)})")

    # Cargas: NAGL (qualidade ~AM1-BCC, roda no Windows) -> AM1-BCC nativo
    # (AmberTools/OpenEye, se houver) -> Gasteiger (stopgap).
    charge_method = None
    try:                                   # 1) NAGL
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        ntkw = NAGLToolkitWrapper()
        model = None
        try:
            from openff.nagl_models import get_models_by_type
            cands = list(get_models_by_type("am1bcc"))
            if cands:
                model = sorted(str(p) for p in cands)[-1]
        except Exception:
            pass
        if model is None:
            for _m in ("openff-gnn-am1bcc-0.1.0-rc.3.pt", "openff-gnn-am1bcc-0.1.0-rc.2.pt",
                       "openff-gnn-am1bcc-0.1.0-rc.1.pt"):
                model = _m
                break
        mname = Path(model).name
        mol_off.assign_partial_charges(mname, toolkit_registry=ntkw)
        charge_method = f"NAGL ({mname})"
        print(f"  [carga] {charge_method} OK  (qualidade ~AM1-BCC)")
    except Exception as e:
        print(f"  [carga] NAGL indisponivel ({str(e)[:110]})")
    if charge_method is None:               # 2) AM1-BCC nativo
        try:
            mol_off.assign_partial_charges("am1bcc")
            charge_method = "AM1-BCC"; print("  [carga] AM1-BCC OK")
        except Exception:
            pass
    if charge_method is None:               # 3) Gasteiger
        from rdkit.Chem import AllChem
        AllChem.ComputeGasteigerCharges(mol_rdkit)
        g = [float(a.GetPropsAsDict().get("_GasteigerCharge", 0.0)) for a in mol_rdkit.GetAtoms()]
        s = sum(g); g = [x - s / len(g) for x in g]
        mol_off.partial_charges = np.array(g) * ou.elementary_charge
        charge_method = "Gasteiger (stopgap)"; print("  [carga] Gasteiger (stopgap)")
    charges = [float(c.m_as(ou.elementary_charge)) for c in mol_off.partial_charges]

    labels, used = None, None
    for sm in ["openff-2.2.1", "openff-2.2.0", "openff-2.1.1", "openff-2.1.0", "openff-2.0.0"]:
        try:
            ff = OFFff(f"{sm}.offxml")
            labels = ff.label_molecules(OFFTop.from_molecules([mol_off]))[0]
            used = sm; print(f"  [SMIRNOFF] labels via {sm}"); break
        except Exception:
            continue
    if labels is None:
        raise RuntimeError("Nenhum forcefield SMIRNOFF disponivel p/ label_molecules")

    # prefixo de tipo derivado do residuo (evita colisao FAD0.. x LIG0.. quando
    # os dois XMLs sao carregados juntos no md_clean). FAD continua "FAD0.." (sem regressao).
    types = [f"{residue_name}{i}" for i in range(n)]
    root = ET.Element("ForceField")
    ET.SubElement(root, "Info").text = (
        f"FAD regen XTAL+SMIRNOFF {used} (charges={charge_method}); {SCRIPT_VERSION}")

    at = ET.SubElement(root, "AtomTypes")
    for i, t in enumerate(types):
        a = mol_off.atoms[i]
        ET.SubElement(at, "Type", name=t, element=a.symbol,
                      mass=str(round(float(a.mass.magnitude), 4)), **{"class": t})

    res = ET.SubElement(ET.SubElement(root, "Residues"), "Residue", name=residue_name)
    for i, (t, nm) in enumerate(zip(types, names)):
        ET.SubElement(res, "Atom", name=nm, type=t, charge=str(round(charges[i], 6)))
    for b in mol_off.bonds:
        ET.SubElement(res, "Bond", atomName1=names[b.atom1_index], atomName2=names[b.atom2_index])

    bf = ET.SubElement(root, "HarmonicBondForce"); seen = set()
    for (i, j), p in labels.get("Bonds", {}).items():
        k = (min(i, j), max(i, j))
        if k in seen: continue
        seen.add(k)
        ET.SubElement(bf, "Bond", type1=types[i], type2=types[j],
                      length=str(round(float(p.length.magnitude) * ANG_NM, 6)),
                      k=str(round(float(p.k.magnitude) * KCAL_KJ / ANG_NM ** 2, 4)))

    af = ET.SubElement(root, "HarmonicAngleForce"); seen = set()
    for (i, j, k), p in labels.get("Angles", {}).items():
        if (i, j, k) in seen: continue
        seen.add((i, j, k))
        ET.SubElement(af, "Angle", type1=types[i], type2=types[j], type3=types[k],
                      angle=str(round(math.radians(float(p.angle.magnitude)), 6)),
                      k=str(round(float(p.k.magnitude) * KCAL_KJ, 4)))

    tf = ET.SubElement(root, "PeriodicTorsionForce"); npr = nim = 0
    for kind, tag in (("ProperTorsions", "Proper"), ("ImproperTorsions", "Improper")):
        for (i, j, k, l), plist in labels.get(kind, {}).items():
            for p in (plist if isinstance(plist, list) else [plist]):
                attrs = {"type1": types[i], "type2": types[j], "type3": types[k], "type4": types[l]}
                for m in range(len(p.periodicity)):
                    div = 3.0 if tag == "Improper" else 1.0
                    attrs[f"periodicity{m+1}"] = str(int(p.periodicity[m]))
                    attrs[f"phase{m+1}"] = str(round(math.radians(float(p.phase[m].magnitude)), 6))
                    attrs[f"k{m+1}"] = str(round(float(p.k[m].magnitude) * KCAL_KJ / div, 4))
                ET.SubElement(tf, tag, **attrs)
                if tag == "Proper": npr += 1
                else: nim += 1

    nf = ET.SubElement(root, "NonbondedForce", coulomb14scale="0.833333", lj14scale="0.5")
    for i, t in enumerate(types):
        vdw = labels.get("vdW", {}).get((i,))
        if vdw is not None:
            eps = float(vdw.epsilon.magnitude) * KCAL_KJ
            sig = 2.0 * float(vdw.rmin_half.magnitude) * ANG_NM / (2.0 ** (1.0 / 6.0))
        else:
            eps, sig = 0.0, 0.35
        ET.SubElement(nf, "Atom", type=t, charge=str(round(charges[i], 6)),
                      sigma=str(round(sig, 6)), epsilon=str(round(eps, 6)))

    print(f"  [SMIRNOFF] bonds={len(labels.get('Bonds', {}))} angles={len(labels.get('Angles', {}))} "
          f"propers={npr} impropers={nim}")
    return '<?xml version="1.0"?>\n' + ET.tostring(root, encoding="unicode")


# ─────────────────────────────────────────────────────────────────────────────
# Gravar o template PDB v3 (86 atomos + CONECT correto)
# ─────────────────────────────────────────────────────────────────────────────
def write_v3_pdb(molH, names):
    conf = molH.GetConformer()
    lines = []
    for i, a in enumerate(molH.GetAtoms()):
        p = conf.GetAtomPosition(i)
        el = a.GetSymbol()
        nm = names[i]
        namecol = f" {nm:<3}" if len(nm) < 4 else nm[:4]
        lines.append(f"HETATM{i+1:5d} {namecol} FAD A   1    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2}")
    # CONECT (1-based)
    nbrs = {i: [] for i in range(molH.GetNumAtoms())}
    for b in molH.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        nbrs[a1].append(a2); nbrs[a2].append(a1)
    for i in range(molH.GetNumAtoms()):
        if nbrs[i]:
            row = "CONECT" + f"{i+1:5d}" + "".join(f"{j+1:5d}" for j in sorted(nbrs[i]))
            lines.append(row)
    lines.append("END")
    OUT_PDB.write_text("\n".join(lines) + "\n")
    print(f"  [pdb] {OUT_PDB.name} gravado ({molH.GetNumAtoms()} atomos, "
          f"{sum(len(v) for v in nbrs.values())//2} ligacoes em CONECT)")


# ─────────────────────────────────────────────────────────────────────────────
def self_validate():
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    print("\n=== AUTO-VALIDACAO (forcefield COMPLETA, igual a MD) ===")
    try:
        # Reproduz o contexto do md_clean: amber14 + tip3p + FAD_v3 + LIG.
        lig = SMOKE / "lig_TcTR_smirnoff.xml"
        ff_files = ["amber14-all.xml", "amber14/tip3p.xml", str(OUT_XML)]
        if lig.exists():
            ff_files.append(str(lig))
        ff = app.ForceField(*ff_files)
        pdb = app.PDBFile(str(OUT_PDB))
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=None, rigidWater=False)
    except Exception as e:
        print(f"  FALHA createSystem (mesmo erro que a MD veria): {str(e)[:200]}")
        return False
    masses = {a.name: system.getParticleMass(a.index).value_in_unit(unit.dalton)
              for a in pdb.topology.atoms()}
    ntor = sum(f.getNumTorsions() for f in system.getForces()
               if isinstance(f, mm.PeriodicTorsionForce))
    print(f"  createSystem OK | torsions={ntor} | N3A={masses.get('N3A',-1):.1f} "
          f"C4A={masses.get('C4A',-1):.1f} PA={masses.get('PA',-1):.1f}")
    integ = mm.LangevinMiddleIntegrator(5 * unit.kelvin, 1 / unit.picosecond,
                                        0.00025 * unit.picoseconds)
    sim = app.Simulation(pdb.topology, system, integ, mm.Platform.getPlatformByName('CPU'))
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=300)
    for _ in range(8):
        sim.step(50)
        e = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        if not np.isfinite(e):
            print("  WARM-UP FAD isolado: NaN"); return False
    ok = ntor > 0 and abs(masses.get('PA', 0) - 31) < 2
    print(f"  WARM-UP FAD isolado (400 passos @5K): OK (E={e:.3e} kJ/mol)")
    print("  RESULTADO:", "PASS" if ok else "ATENCAO (revise torsions/massas)")
    return ok


def main():
    print(f"[regen_fad_xml.py {SCRIPT_VERSION}]  (fonte: cristal -> CCD -> SMIRNOFF)")
    print(f"XTAL  : {XTAL_PDB}")
    print(f"Saidas: {OUT_PDB.name} + {OUT_XML.name}\n")
    if not XTAL_PDB.exists():
        sys.exit(f"FALTA {XTAL_PDB}")

    smiles = load_ccd_smiles()
    print(f"  [CCD] SMILES: {smiles[:70]}{'...' if smiles and len(smiles) > 70 else ''}")
    molH, names = build_fad_from_xtal(smiles)
    write_v3_pdb(molH, names)
    xml = build_xml(molH, names, residue_name="FAD")
    OUT_XML.write_text(xml, encoding="utf-8")
    print(f"\nGravado: {OUT_XML.name}")

    ok = self_validate()

    # Sincroniza os arquivos ORIGINAIS (que o md_clean usa por default) com o v3,
    # SO se a auto-validacao passou (espelha regen_ligand_xml.py). Assim nao
    # clobbera um original bom quando a regeneracao falha.
    # Evita o caso topologia-v3 (correta) + force-field-antigo (embaralhado) ->
    # "No template found ... bonds are different". Faz backup do original 1x.
    if ok:
        import shutil
        for _v3, _orig in [(OUT_XML, SMOKE / "cofactor_FAD_TcTR.xml"),
                           (OUT_PDB, SMOKE / "cofactor_FAD_TcTR.pdb")]:
            _bak = _orig.with_suffix(_orig.suffix + ".broken_bak")
            if _orig.exists() and not _bak.exists():
                shutil.copy2(_orig, _bak)
            shutil.copy2(_v3, _orig)
            print(f"  [sync] {_orig.name} <- {_v3.name}")
    else:
        print("  [sync] PULADO (validacao falhou; nao sobrescrevo os originais).")

    print("\n" + "=" * 64)
    print("PROXIMO PASSO: rodar a MD (o md_simulation.py ja prefere os arquivos *_v3.*):")
    print('  python md_simulation.py --fresh --ns 1 --tag prod_v3 --platform CUDA')
    print("Confira no log: 'usando FAD XML regenerado' e SEM 'BOND ANOMALA'.")
    print("NOTA: para producao/publicacao prefira cargas AM1-BCC (AmberTools sqm).")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
