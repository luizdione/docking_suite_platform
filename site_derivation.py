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
    import docking_v2
    from Bio.PDB import MMCIFParser, PDBParser, Superimposer
    from Bio.Align import PairwiseAligner, substitution_matrices

    resname = parse_resname(reference_ligand_resname)

    # alvo (AlphaFold .cif) -- 1a cadeia do 1o modelo
    tgt = MMCIFParser(QUIET=True).get_structure("tgt", str(target_cif))
    tgt_chain = next(iter(next(iter(tgt))))
    tgt_cas, tgt_seq = _chain_ca_seq(tgt_chain)
    if len(tgt_cas) < min_ca_matched:
        raise RuntimeError(f"alvo {target_cif} tem poucos CAs ({len(tgt_cas)})")

    # referencia (co-cristal do RCSB)
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

    # copia do ligante mais proxima da cadeia proteica alinhada
    def _cen(res):
        return np.mean([a.get_coord() for a in res.get_atoms()], axis=0)
    lig_res, _cid = min(copies, key=lambda rc: np.linalg.norm(_cen(rc[0]) - ref_ca_centroid))
    lig_coords = np.array([a.get_coord() for a in lig_res.get_atoms()])

    # alinhamento de sequencia -> pares de CA 1:1
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aln = aligner.align(tgt_seq, ref_seq)[0]
    fixed, moving = [], []  # fixed=alvo, moving=referencia
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
    rmsd = float(sup.rms)
    if rmsd > rmsd_max_A:
        raise RuntimeError(
            f"superposicao ruim contra {reference_pdb}: RMSD {rmsd:.2f} A > "
            f"{rmsd_max_A} A sobre {n} CAs -- centro NAO confiavel (fail-loud)")

    center_crystal = lig_coords.mean(axis=0)                 # convencao A
    center_target = (np.dot(lig_coords, rot) + tran).mean(axis=0)  # convencao B

    extent = lig_coords.max(axis=0) - lig_coords.min(axis=0)
    box = np.maximum(extent + 2 * padding_A, min_box_A)

    return {
        "reference_pdb": reference_pdb,
        "reference_ligand_resname": resname,
        "n_ligand_atoms": int(len(lig_coords)),
        "n_ca_matched": int(n),
        "superposition_rmsd_A": round(rmsd, 3),
        "box_center_A_target_frame": [round(float(x), 3) for x in center_target],
        "box_center_A_crystal_frame": [round(float(x), 3) for x in center_crystal],
        "box_size_A": [round(float(x), 1) for x in box],
    }
