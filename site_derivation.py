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


def _extract_reference_ligand(structure, reference_ligand_resname, ref_chain_centroid):
    """Extrai o ligante de referencia. Dois modos:
      * HETATM de 1-3 letras (ex.: 'CY9', 'LAR'): copia mais proxima do centroide
        da cadeia proteica de referencia (logica intraprotomerica).
      * cadeia peptidica (peptideo ciclico modelado como polimero, ex.: faloidina em
        7BTI: 'cadeia polimerica X (7 res...)'): junta os atomos da cadeia-copia mais
        proxima do centroide. Antes isto caia em parse_resname()->'CADEIA' e falhava
        com RuntimeError; agora vira um ligante real (e o gate de imersao decide).
    Devolve (coords, elementos, modo, rotulo). Falha ALTO se nao achar -- nao inventar."""
    field = str(reference_ligand_resname)
    m = re.search(r"cadeia\s+polimerica\s+([A-Za-z0-9])", field, re.IGNORECASE)
    if m:
        copies_ids = re.findall(r"'([A-Za-z0-9])'", field) or [m.group(1)]
        model = next(iter(structure))
        chain_ids = [c.id for c in model]
        best, best_cid, best_d = None, None, 1e18
        for cid in copies_ids:
            if cid not in chain_ids:
                continue
            atoms = [a for res in model[cid] for a in res.get_atoms()]
            if not atoms:
                continue
            cen = np.mean([a.get_coord() for a in atoms], axis=0)
            d = float(np.linalg.norm(cen - ref_chain_centroid))
            if d < best_d:
                best, best_cid, best_d = atoms, cid, d
        if best is None:
            raise RuntimeError(
                f"cadeia peptidica {copies_ids} (de {field!r}) nao encontrada -- "
                f"nao inventar centro")
        coords = np.array([a.get_coord() for a in best])
        elems = [a.element.strip().capitalize() for a in best]
        return coords, elems, "peptide_chain", f"PEP:{best_cid}"

    resname = parse_resname(field)
    copies = _find_ligand_copies(structure, resname)
    if not copies:
        raise RuntimeError(
            f"ligante '{resname}' (de {field!r}) nao encontrado em referencia -- "
            f"nao inventar centro")

    def _cen(res):
        return np.mean([a.get_coord() for a in res.get_atoms()], axis=0)
    lig_res, _cid = min(copies, key=lambda rc: np.linalg.norm(_cen(rc[0]) - ref_chain_centroid))
    atoms = list(lig_res.get_atoms())
    coords = np.array([a.get_coord() for a in atoms])
    elems = [a.element.strip().capitalize() for a in atoms]
    return coords, elems, "hetatm", resname


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
    sequencia + Kabsch, e localiza o ligante de referencia (HETATM ou cadeia
    peptidica). NAO aplica limiar de RMSD (o chamador decide). Devolve rot/tran,
    RMSD, n_ca, o ligante (coords do cristal + elementos + modo), a estrutura do
    alvo (p/ calcular imersao) e os PARES de CA alinhados fixed_ca/moving_ca (p/ o
    chamador tentar superposicao LOCAL do bolso quando a global for ruim)."""
    import docking_v2
    from Bio.PDB import MMCIFParser, PDBParser, Superimposer
    from Bio.Align import PairwiseAligner, substitution_matrices

    tgt = MMCIFParser(QUIET=True).get_structure("tgt", str(target_cif))
    tgt_chain = next(iter(next(iter(tgt))))
    tgt_cas, tgt_seq = _chain_ca_seq(tgt_chain)
    if len(tgt_cas) < min_ca_matched:
        raise RuntimeError(f"alvo {target_cif} tem poucos CAs ({len(tgt_cas)})")

    pdb_path, msg = docking_v2.download_pdb(reference_pdb, download_dir)
    if pdb_path is None:
        raise RuntimeError(f"nao baixei {reference_pdb}: {msg}")
    ref = PDBParser(QUIET=True).get_structure("ref", str(pdb_path))

    ref_chain = _biggest_protein_chain(ref)
    ref_cas, ref_seq = _chain_ca_seq(ref_chain)
    ref_ca_centroid = np.mean([a.get_coord() for a in ref_cas], axis=0)

    lig_coords, lig_elems, lig_mode, lig_label = _extract_reference_ligand(
        ref, reference_ligand_resname, ref_ca_centroid)

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
        "n_ca_matched": int(n), "resname": lig_label, "lig_mode": lig_mode,
        "lig_coords": lig_coords, "lig_elements": lig_elems,
        "target_structure": tgt, "fixed_ca": fixed, "moving_ca": moving,
    }


def geometric_validation(target_cif, reference_pdb, reference_ligand_resname,
                         download_dir, contact_A=4.5, min_frac_contact=0.8,
                         max_mean_A=4.0, rmsd_max_A=3.0, min_ca_matched=30,
                         local_rmsd_max_A=2.5, pocket_r_A=12.0, min_local_ca=15):
    """Valida GEOMETRICAMENTE (sem docking, sem GPU) um sitio derivado: a pose de
    referencia transformada p/ o frame do alvo cai num BOLSO real do modelo?

    Redock rigido nao serve p/ estes sitios (modelos AlphaFold apo -- bolso pode
    estar fechado; o Vina otimiza energia e foge do bolso). O gate honesto e este:
    fracao de atomos pesados do ligante com um atomo de proteina < contact_A, e a
    distancia media ao receptor. Discrimina bolso presente x ausente (ex.: ligante
    de FILAMENTO num modelo MONOMERICO da imersao fraca -- biologicamente correto).

    RESGATE por superposicao LOCAL (paralogos): quando a superposicao GLOBAL passa de
    rmsd_max_A (dobra global divergente -- ex.: TcAct2-4 com insercoes longe do bolso),
    tenta um Kabsch so nos CAs que forram o bolso (ref-CA a <=pocket_r_A do ligante).
    Se o ajuste LOCAL for bom (<=local_rmsd_max_A sobre >=min_local_ca CAs) E a pose
    continuar IMERSA nos atomos REAIS do alvo -> veredito 'center_ok_local'. A imersao
    e o backstop honesto: ligante de filamento consegue ajuste local otimo mas NAO
    imerge num monomero, logo NAO e resgatado (permanece weak_immersion)."""
    from scipy.spatial import cKDTree
    from Bio.PDB import Superimposer
    g = _superpose_target_reference(target_cif, reference_pdb, reference_ligand_resname,
                                    download_dir, min_ca_matched=min_ca_matched)
    rmsd = g["superposition_rmsd_A"]
    out = {"reference_pdb": reference_pdb, "reference_ligand_resname": g["resname"],
           "superposition_rmsd_A": rmsd, "n_ca_matched": g["n_ca_matched"],
           "ligand_mode": g.get("lig_mode", "hetatm")}

    coords, elems = g["lig_coords"], g["lig_elements"]
    heavy = np.array([c for c, e in zip(coords, elems) if e.upper() != "H"])
    rec = np.array([a.get_coord() for a in g["target_structure"].get_atoms()
                    if a.element.strip().upper() != "H"])
    rec_tree = cKDTree(rec)

    def _immersion(rot, tran):
        expected = np.dot(heavy, rot) + tran                    # pose no frame do alvo
        d, _ = rec_tree.query(expected)
        box = np.maximum((expected.max(axis=0) - expected.min(axis=0)) + 16.0, 22.5)
        return {
            "immersion_mean_A": round(float(d.mean()), 2),
            "immersion_min_A": round(float(d.min()), 2),
            "frac_contact": round(float((d < contact_A).mean()), 2),
            "box_center_A_target_frame": [round(float(x), 3) for x in expected.mean(axis=0)],
            "box_size_A": [round(float(x), 1) for x in box],
        }

    def _immersed(imm):
        return imm["frac_contact"] >= min_frac_contact and imm["immersion_mean_A"] < max_mean_A

    if rmsd <= rmsd_max_A:                                       # superposicao global confiavel
        imm = _immersion(g["rot"], g["tran"])
        out.update(imm)
        out["superposition_mode"] = "global"
        out["verdict"] = "center_ok" if _immersed(imm) else "weak_immersion"
        return out

    # --- global ruim: RESGATE por superposicao local do bolso ---
    moving_xyz = np.array([a.get_coord() for a in g["moving_ca"]])   # ref-CAs (frame do cristal)
    dmin, _ = cKDTree(heavy).query(moving_xyz)
    keep = np.where(dmin <= pocket_r_A)[0]
    out["n_ca_local"] = int(len(keep))
    if len(keep) < min_local_ca:
        out["superposition_mode"] = "local_pocket"
        out["verdict"] = "blocked_superposition"
        return out
    sup = Superimposer()
    sup.set_atoms([g["fixed_ca"][i] for i in keep], [g["moving_ca"][i] for i in keep])
    l_rmsd = round(float(sup.rms), 3)
    l_rot, l_tran = sup.rotran
    imm = _immersion(l_rot, l_tran)
    out.update(imm)
    out["local_rmsd_A"] = l_rmsd
    out["superposition_mode"] = "local_pocket"
    if l_rmsd <= local_rmsd_max_A and _immersed(imm):
        out["verdict"] = "center_ok_local"
    elif l_rmsd <= local_rmsd_max_A:
        out["verdict"] = "weak_immersion"        # bolso alinhavel mas ligante nao imerso
    else:
        out["verdict"] = "blocked_superposition"  # nem local ajusta -> divergencia real
    return out


def validate_embedded_ligand(receptor_file, ligand_resname, contact_A=4.5,
                             min_frac_contact=0.8, max_mean_A=4.0, padding_A=8.0,
                             min_box_A=22.5):
    """Valida um sitio cujo RECEPTOR JA CONTEM o ligante na pose nativa -- ex.: o
    pentamero de F-actina com jasplaquinolida (9UE) embutida (threading sobre 7Q8C no
    frame do 7PM3). Diferente de geometric_validation, aqui NAO ha superposicao: o
    ligante ja esta no referencial do receptor. Multi-cadeia por construcao (usa TODAS
    as cadeias do receptor), o que resolve o sitio INTER-PROTOMERICO que a superposicao
    de cadeia unica nao alcanca.

    Se ha varias copias do ligante (o 5-mer tem 3 copias de 9UE, uma por juncao),
    escolhe a MAIS IMERSA -- a copia CENTRAL, totalmente enterrada; nao a de borda,
    parcialmente exposta (foi a borda ~36% exposta que fez o TRIMERO ser rejeitado a
    favor do pentamero). Devolve o centro/box dessa copia.

    Metrica identica a geometric_validation: fracao de atomos pesados do ligante a
    <contact_A de um atomo do receptor, e distancia media. Veredito center_ok/
    weak_immersion. Fail-loud (RuntimeError) se o ligante nao existe no arquivo --
    nunca inventar centro."""
    from scipy.spatial import cKDTree
    from Bio.PDB import PDBParser, MMCIFParser
    p = str(receptor_file)
    parser = MMCIFParser(QUIET=True) if p.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    struct = parser.get_structure("rec", p)

    lig_copies = {}           # (chain_id, resseq) -> [heavy atoms]
    rec_atoms = []            # atomos pesados do receptor (tudo menos o proprio ligante)
    target = ligand_resname.strip().upper()
    for model in struct:
        for chain in model:
            for res in chain:
                is_lig = res.get_resname().strip().upper() == target
                heavy = [a for a in res.get_atoms() if a.element.strip().upper() != "H"]
                if is_lig:
                    if heavy:
                        lig_copies[(chain.id, res.id[1])] = heavy
                else:
                    rec_atoms.extend(a.get_coord() for a in heavy)
        break                 # so o 1o modelo
    if not lig_copies:
        raise RuntimeError(
            f"ligante '{ligand_resname}' nao encontrado em {receptor_file} -- nao inventar centro")
    if not rec_atoms:
        raise RuntimeError(f"receptor {receptor_file} sem atomos (fora o ligante)")

    tree = cKDTree(np.array(rec_atoms))
    best = None
    for (cid, rseq), atoms in lig_copies.items():
        coords = np.array([a.get_coord() for a in atoms])
        d, _ = tree.query(coords)
        frac = float((d < contact_A).mean())
        mean_d = float(d.mean())
        ext = coords.max(axis=0) - coords.min(axis=0)
        box = np.maximum(ext + 2 * padding_A, min_box_A)
        rec = {
            "copy_chain": cid, "copy_resseq": int(rseq), "n_ligand_heavy": len(atoms),
            "frac_contact": round(frac, 2), "immersion_mean_A": round(mean_d, 2),
            "immersion_min_A": round(float(d.min()), 2),
            "box_center_A": [round(float(x), 3) for x in coords.mean(axis=0)],
            "box_size_A": [round(float(x), 1) for x in box],
            "_score": (frac, -mean_d),
        }
        if best is None or rec["_score"] > best["_score"]:
            best = rec
    best.pop("_score")
    best["n_ligand_copies"] = len(lig_copies)
    best["reference_ligand_resname"] = target
    best["box_center_frame"] = "receptor_native"
    best["verdict"] = ("center_ok" if (best["frac_contact"] >= min_frac_contact
                       and best["immersion_mean_A"] < max_mean_A) else "weak_immersion")
    return best
