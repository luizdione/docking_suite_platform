# -*- coding: utf-8 -*-
"""site_derivation.py -- deriva o centro do sitio de ligacao por SUPERPOSICAO
estrutural contra um complexo de referencia (co-cristal).

Motivacao (SUB04): alvos AlphaFold puros (actinas de T. cruzi) nao tem ligante
co-cristalizado proprio. O sitio de uma droga do citoesqueleto e definido por um
complexo de referencia experimental (ex.: latrunculina em 1ESV, citocalasina em
3EKS). Este modulo alinha a sequencia do alvo com a da referencia, superpoe os
CAs (Kabsch/SVD via Bio.PDB.Superimposer) e transporta o centroide do ligante de
referencia para derivar o centro da caixa de docking.

DOIS REFERENCIAIS (importante -- nao confundir; foi o que levou 88 A de erro no
primeiro prototipo):
  * convencao A ("crystal_frame"): centro = centroide do ligante nas coordenadas
    BRUTAS do co-cristal. So faz sentido com o receptor MOVIDO para o frame do
    cristal (prepared_receptor_local_file). E a convencao do unico par ja validado
    (TcAct1 x latrunculina-A, explicit_center).
  * convencao B ("target_frame"): centro = centroide do ligante TRANSFORMADO para
    o referencial do proprio alvo (AlphaFold). Usa o receptor AlphaFold INTOCADO.
    Mais simples (nao gera receptor derivado), recomendado para generalizar.

Este modulo COMPUTA e RETORNA os dois; NAO escreve no registry e NAO decide a
convencao -- quem persiste (tools/derive_sites.py) decide, e so depois de um
redock gate confirmar. Falha ALTO (RuntimeError) quando nao da para confiar:
ligante nao encontrado, poucos CAs pareados, ou RMSD de superposicao alto.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def three_to_one(resname):
    """3->1 letra, robusto entre versoes de Biopython (three_to_one foi removido)."""
    r = resname.strip()
    for mod in ("Bio.Data.PDBData", "Bio.Data.IUPACData"):
        try:
            d = __import__(mod, fromlist=["protein_letters_3to1"]).protein_letters_3to1
            if r.upper() in d:
                return d[r.upper()]
            if r.capitalize() in d:
                return d[r.capitalize()]
        except Exception:
            pass
    return "X"


def parse_resname(raw):
    """Extrai o codigo do ligante de um campo possivelmente 'sujo'.

    O registry guarda coisas como 'HET CY9 copia na cadeia A (94 % enterrada...)'.
    O codigo real e o token de 1-3 caracteres logo apos 'HET', ou o proprio campo
    se ja vier limpo (ex.: 'LAR', 'LAB').
    """
    s = str(raw).strip()
    m = re.match(r"^HET\s+([A-Za-z0-9]{1,3})\b", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    tok = s.split()[0] if s.split() else s
    return re.sub(r"[^A-Za-z0-9]", "", tok).upper()


def _chain_ca_seq(chain):
    """(lista de CA atoms, sequencia 1-letra) dos residuos padrao com CA."""
    cas, seq = [], []
    for res in chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        cas.append(res["CA"])
        seq.append(three_to_one(res.get_resname()))
    return cas, "".join(seq)


def _find_ligand_copies(structure, resname3):
    """Todas as copias (residuos HETATM) cujo resname bate: [(residuo, chain_id)]."""
    out = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip().upper() == resname3.upper():
                    out.append((res, chain.id))
        break  # so o 1o modelo
    return out


def _biggest_protein_chain(structure):
    """Cadeia com mais CAs (representa a proteina principal da referencia)."""
    best, best_n = None, -1
    for model in structure:
        for chain in model:
            cas, _ = _chain_ca_seq(chain)
            if len(cas) > best_n:
                best, best_n = chain, len(cas)
        break
    return best


def derive_site_center(target_cif, reference_pdb, reference_ligand_resname,
                       download_dir, rmsd_max_A=3.0, min_ca_matched=30,
                       padding_A=8.0, min_box_A=22.5):
    """Deriva o centro do sitio por superposicao alvo<-referencia.

    Retorna dict com box_center_A_target_frame (convencao B),
    box_center_A_crystal_frame (convencao A), box_size_A, superposition_rmsd_A,
    n_ca_matched, reference_pdb, reference_ligand_resname (limpo). Levanta
    RuntimeError se algo impede confiar no resultado.
    """
    g = _superpose_target_reference(target_cif, reference_pdb, reference_ligand_resname,
                                    download_dir, min_ca_matched=min_ca_matched)
    n, rmsd = g["n_ca_matched"], g["superposition_rmsd_A"]
    if rmsd > rmsd_max_A:
        raise RuntimeError(
            f"superposicao ruim contra {reference_pdb}: RMSD {rmsd:.2f} A > "
            f"{rmsd_max_A} A sobre {n} CAs -- centro NAO confiavel (fail-loud)")

    lig_coords, rot, tran = g["lig_coords"], g["rot"], g["tran"]
    center_crystal = lig_coords.mean(axis=0)                       # convencao A
    center_target = (np.dot(lig_coords, rot) + tran).mean(axis=0)  # convencao B
    extent = lig_coords.max(axis=0) - lig_coords.min(axis=0)
    box = np.maximum(extent + 2 * padding_A, min_box_A)

    return {
        "reference_pdb": reference_pdb,
        "reference_ligand_resname": g["resname"],
        "n_ligand_atoms": int(len(lig_coords)),
        "n_ca_matched": int(n),
        "superposition_rmsd_A": round(rmsd, 3),
        "box_center_A_target_frame": [round(float(x), 3) for x in center_target],
        "box_center_A_crystal_frame": [round(float(x), 3) for x in center_crystal],
        "box_size_A": [round(float(x), 1) for x in box],
    }


def _superpose_target_reference(target_cif, reference_pdb, reference_ligand_resname,
                                download_dir, min_ca_matched=30):
    """Superpoe o alvo (AlphaFold) na referencia (co-cristal) via alinhamento de
    sequencia + Kabsch, e localiza o ligante de referencia. NAO aplica limiar de
    RMSD (o chamador decide). Devolve rot/tran, RMSD, n_ca, o ligante (coords do
    cristal + elementos) e a estrutura do alvo (p/ calcular imersao)."""
    import docking_v2
    from Bio.PDB import MMCIFParser, PDBParser, Superimposer
    from Bio.Align import PairwiseAligner, substitution_matrices

    resname = parse_resname(reference_ligand_resname)
    tgt = MMCIFParser(QUIET=True).get_structure("tgt", str(target_cif))
    tgt_chain = next(iter(next(iter(tgt))))
    tgt_cas, tgt_seq = _chain_ca_seq(tgt_chain)
    if len(tgt_cas) < min_ca_matched:
        raise RuntimeError(f"alvo {target_cif} tem poucos CAs ({len(tgt_cas)})")

    pdb_path, msg = docking_v2.download_pdb(reference_pdb, download_dir)
    if pdb_path is None:
        raise RuntimeError(f"nao baixei {reference_pdb}: {msg}")
    ref = PDBParser(QUIET=True).get_structure("ref", str(pdb_path))

    copies = _find_ligand_copies(ref, resname)
    if not copies:
        raise RuntimeError(
            f"ligante '{resname}' (de {reference_ligand_resname!r}) nao encontrado "
            f"em {reference_pdb} -- nao inventar centro")

    ref_chain = _biggest_protein_chain(ref)
    ref_cas, ref_seq = _chain_ca_seq(ref_chain)
    ref_ca_centroid = np.mean([a.get_coord() for a in ref_cas], axis=0)

    def _cen(res):
        return np.mean([a.get_coord() for a in res.get_atoms()], axis=0)
    lig_res, _cid = min(copies, key=lambda rc: np.linalg.norm(_cen(rc[0]) - ref_ca_centroid))
    lig_atoms = list(lig_res.get_atoms())
    lig_coords = np.array([a.get_coord() for a in lig_atoms])
    lig_elems = [a.element.strip().capitalize() for a in lig_atoms]

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aln = aligner.align(tgt_seq, ref_seq)[0]
    fixed, moving = [], []
    for (t0, t1), (r0, r1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(t1 - t0):
            fixed.append(tgt_cas[t0 + k])
            moving.append(ref_cas[r0 + k])
    n = len(fixed)
    if n < min_ca_matched:
        raise RuntimeError(f"poucos CAs pareados ({n}) entre alvo e {reference_pdb}")

    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    rot, tran = sup.rotran
    return {
        "rot": rot, "tran": tran, "superposition_rmsd_A": round(float(sup.rms), 3),
        "n_ca_matched": int(n), "resname": resname,
        "lig_coords": lig_coords, "lig_elements": lig_elems,
        "target_structure": tgt,
    }


def geometric_validation(target_cif, reference_pdb, reference_ligand_resname,
                         download_dir, contact_A=4.5, min_frac_contact=0.8,
                         max_mean_A=4.0, rmsd_max_A=3.0, min_ca_matched=30):
    """Valida GEOMETRICAMENTE (sem docking, sem GPU) um sitio derivado: a pose de
    referencia transformada p/ o frame do alvo cai num BOLSO real do modelo?

    Redock rigido nao serve p/ estes sitios (modelos AlphaFold apo -- bolso pode
    estar fechado; o Vina otimiza energia e foge do bolso). O gate honesto e este:
    fracao de atomos pesados do ligante com um atomo de proteina < contact_A, e a
    distancia media ao receptor. Discrimina bolso presente x ausente (ex.: ligante
    de FILAMENTO num modelo MONOMERICO da imersao fraca -- biologicamente correto)."""
    from scipy.spatial import cKDTree
    g = _superpose_target_reference(target_cif, reference_pdb, reference_ligand_resname,
                                    download_dir, min_ca_matched=min_ca_matched)
    rmsd = g["superposition_rmsd_A"]
    out = {"reference_pdb": reference_pdb, "reference_ligand_resname": g["resname"],
           "superposition_rmsd_A": rmsd, "n_ca_matched": g["n_ca_matched"]}
    if rmsd > rmsd_max_A:
        out["verdict"] = "blocked_superposition"
        return out

    coords, elems = g["lig_coords"], g["lig_elements"]
    heavy = np.array([c for c, e in zip(coords, elems) if e.upper() != "H"])
    expected = np.dot(heavy, g["rot"]) + g["tran"]              # pose no frame do alvo
    rec = np.array([a.get_coord() for a in g["target_structure"].get_atoms()
                    if a.element.strip().upper() != "H"])
    d, _ = cKDTree(rec).query(expected)
    frac = float((d < contact_A).mean())
    mean_d = float(d.mean())
    box = np.maximum((expected.max(axis=0) - expected.min(axis=0)) + 16.0, 22.5)
    out.update({
        "immersion_mean_A": round(mean_d, 2),
        "immersion_min_A": round(float(d.min()), 2),
        "frac_contact": round(frac, 2),
        "box_center_A_target_frame": [round(float(x), 3) for x in expected.mean(axis=0)],
        "box_size_A": [round(float(x), 1) for x in box],
        "verdict": "center_ok" if (frac >= min_frac_contact and mean_d < max_mean_A)
                   else "weak_immersion",
    })
    return out
