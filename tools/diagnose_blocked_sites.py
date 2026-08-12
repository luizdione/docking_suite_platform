# -*- coding: utf-8 -*-
"""tools/diagnose_blocked_sites.py -- DIAGNOSTICO read-only dos sitios SUB04 que o
gate geometrico global (tools/validate_derived_sites.py) bloqueia. NAO escreve no
registry, NAO altera site_derivation.py, NAO usa GPU/docking. So mede e relata.

Investiga duas hipoteses de por que 13 dos 24 sitios bloqueiam:

  (B) blocked_superposition (9 sitios, TcAct2/3/4): a superposicao Kabsch e GLOBAL
      (~375 CAs). Se um paralogo diverge num subdominio distante do bolso, o RMSD
      global estoura 3 A mesmo com o bolso conservado. Hipotese: superposicao LOCAL
      (so os CAs que forram o bolso, <=pocket_R do ligante de referencia) da um RMSD
      menor E -- o teste honesto -- a pose transportada continua imersa nos atomos
      REAIS do alvo. Imersao nao se "trapaceia" escolhendo poucos pontos de ajuste:
      ou o ligante cai dentro do bolso do alvo, ou nao cai.

  (A) blocked:RuntimeError (4 sitios, faloidina/7BTI): NAO e geometria, e parsing --
      o ligante de referencia e um heptapeptideo ciclico modelado como CADEIA
      polimerica (X/Y/Z), nao um HETATM de 3 letras. parse_resname() pega "CADEIA"
      e _find_ligand_copies() nao acha nada. Aqui extraimos a cadeia peptidica e
      rodamos a MESMA imersao -- esperado: imersao FRACA (faloidina liga FILAMENTO,
      igual jasplaquinolida), i.e. o veredito honesto e weak_immersion, nao um crash.

CHECK DE HONESTIDADE: a superposicao local NAO pode resgatar ligante de filamento
(jasplaquinolida/faloidina) -- o bolso deles exige protomeros vizinhos ausentes no
monomero. Rodamos sobre os 24 sitios: se o metodo local "aprovasse" um filamento,
estaria quebrado.

Uso:
    python tools/diagnose_blocked_sites.py --campaign SUB04
    python tools/diagnose_blocked_sites.py --campaign SUB04 --pocket-r 12 --json out.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry            # noqa: E402
import site_derivation as sd  # noqa: E402  (reusa parse_resname/_chain_ca_seq/etc.)

DL_DIR = ROOT / ".pdb_cache"

# gate local (espelha o global; imersao e o backstop honesto)
CONTACT_A = 4.5
MIN_FRAC = 0.8
MAX_MEAN_A = 4.0
LOCAL_RMSD_MAX_A = 2.5     # ajuste local so vale se for de fato bom
MIN_LOCAL_CA = 15          # e com pontos suficientes


def _iter_derived(campaign_id):
    for target in registry.list_targets(campaign_id):
        tid = target.get("id")
        cif = registry.resolve_local_file((target.get("structure") or {}).get("local_file"))
        for sname, site in (target.get("sites") or {}).items():
            if site.get("definition") == "derived_from_reference_complex":
                yield tid, sname, site, cif


def _aligned_ca_pairs(tgt_cas, tgt_seq, ref_cas, ref_seq):
    """Mesmo alinhamento BLOSUM62 do motor; devolve pares (tgt_CA, ref_CA)."""
    from Bio.Align import PairwiseAligner, substitution_matrices
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
    return fixed, moving


def _extract_reference_ligand(ref, ref_field, ref_chain_centroid):
    """Devolve (coords, elems, mode). Modo HETATM (resname 3 letras) OU modo CADEIA
    polimerica (peptideo: faloidina). Escolhe a copia mais proxima do centroide da
    cadeia proteica de referencia (mesma logica intraprotomerica do motor)."""
    m = re.search(r"cadeia\s+polimerica\s+([A-Za-z0-9])", str(ref_field), re.IGNORECASE)
    if m:
        # peptideo ciclico modelado como cadeia -> junta os atomos da cadeia-copia
        copies_ids = re.findall(r"'([A-Za-z0-9])'", str(ref_field)) or [m.group(1)]
        model = next(iter(ref))
        best, best_d = None, 1e9
        for cid in copies_ids:
            if cid not in [c.id for c in model]:
                continue
            atoms = [a for res in model[cid] for a in res.get_atoms()]
            if not atoms:
                continue
            cen = np.mean([a.get_coord() for a in atoms], axis=0)
            d = float(np.linalg.norm(cen - ref_chain_centroid))
            if d < best_d:
                best, best_d = atoms, d
        if best is None:
            raise RuntimeError(f"cadeia peptidica {copies_ids} nao encontrada")
        coords = np.array([a.get_coord() for a in best])
        elems = [a.element.strip().capitalize() for a in best]
        return coords, elems, "peptide_chain"

    resname = sd.parse_resname(ref_field)
    copies = sd._find_ligand_copies(ref, resname)
    if not copies:
        raise RuntimeError(f"ligante HETATM '{resname}' nao encontrado")

    def _cen(res):
        return np.mean([a.get_coord() for a in res.get_atoms()], axis=0)
    lig_res, _cid = min(copies, key=lambda rc: np.linalg.norm(_cen(rc[0]) - ref_chain_centroid))
    atoms = list(lig_res.get_atoms())
    coords = np.array([a.get_coord() for a in atoms])
    elems = [a.element.strip().capitalize() for a in atoms]
    return coords, elems, "hetatm"


def _kabsch(fixed_atoms, moving_atoms):
    """RMSD + (rot, tran) que leva moving->fixed. Espelha Bio.PDB.Superimposer."""
    from Bio.PDB import Superimposer
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    rot, tran = sup.rotran
    return float(sup.rms), rot, tran


def _immersion(lig_coords, lig_elems, rot, tran, tgt_structure):
    """Transforma o ligante p/ o frame do alvo e mede imersao nos atomos REAIS."""
    from scipy.spatial import cKDTree
    heavy = np.array([c for c, e in zip(lig_coords, lig_elems) if e.upper() != "H"])
    expected = np.dot(heavy, rot) + tran
    rec = np.array([a.get_coord() for a in tgt_structure.get_atoms()
                    if a.element.strip().upper() != "H"])
    d, _ = cKDTree(rec).query(expected)
    return {
        "frac_contact": round(float((d < CONTACT_A).mean()), 2),
        "immersion_mean_A": round(float(d.mean()), 2),
        "immersion_min_A": round(float(d.min()), 2),
        "box_center_A_target_frame": [round(float(x), 3) for x in expected.mean(axis=0)],
    }


def diagnose_site(target_cif, reference_pdb, ref_ligand_field, pocket_r):
    """Global vs local: RMSD e imersao em ambos. Read-only."""
    import docking_v2
    from Bio.PDB import MMCIFParser, PDBParser

    tgt = MMCIFParser(QUIET=True).get_structure("tgt", str(target_cif))
    tgt_chain = next(iter(next(iter(tgt))))
    tgt_cas, tgt_seq = sd._chain_ca_seq(tgt_chain)

    pdb_path, msg = docking_v2.download_pdb(reference_pdb, DL_DIR)
    if pdb_path is None:
        raise RuntimeError(f"nao baixei {reference_pdb}: {msg}")
    ref = PDBParser(QUIET=True).get_structure("ref", str(pdb_path))
    ref_chain = sd._biggest_protein_chain(ref)
    ref_cas, ref_seq = sd._chain_ca_seq(ref_chain)
    ref_centroid = np.mean([a.get_coord() for a in ref_cas], axis=0)

    lig_coords, lig_elems, mode = _extract_reference_ligand(ref, ref_ligand_field, ref_centroid)
    lig_heavy = np.array([c for c, e in zip(lig_coords, lig_elems) if e.upper() != "H"])

    fixed, moving = _aligned_ca_pairs(tgt_cas, tgt_seq, ref_cas, ref_seq)
    if len(fixed) < 30:
        raise RuntimeError(f"poucos CAs pareados ({len(fixed)})")

    # --- global ---
    g_rmsd, g_rot, g_tran = _kabsch(fixed, moving)
    g_imm = _immersion(lig_coords, lig_elems, g_rot, g_tran, tgt)

    # --- local: so os pares cujo ref-CA forra o bolso (<= pocket_r do ligante) ---
    ref_ca_xyz = np.array([a.get_coord() for a in moving])
    from scipy.spatial import cKDTree
    dmin, _ = cKDTree(lig_heavy).query(ref_ca_xyz)  # dist de cada ref-CA ao ligante
    keep = np.where(dmin <= pocket_r)[0]
    n_local = int(len(keep))
    out = {
        "reference_pdb": reference_pdb, "ligand_mode": mode,
        "n_ligand_heavy": int(len(lig_heavy)), "n_ca_global": len(fixed),
        "global_rmsd_A": round(g_rmsd, 3),
        "global_frac_contact": g_imm["frac_contact"],
        "global_immersion_mean_A": g_imm["immersion_mean_A"],
        "n_ca_local": n_local, "pocket_r_A": pocket_r,
    }
    if n_local >= MIN_LOCAL_CA:
        fx = [fixed[i] for i in keep]
        mv = [moving[i] for i in keep]
        l_rmsd, l_rot, l_tran = _kabsch(fx, mv)
        l_imm = _immersion(lig_coords, lig_elems, l_rot, l_tran, tgt)
        out.update({
            "local_rmsd_A": round(l_rmsd, 3),
            "local_frac_contact": l_imm["frac_contact"],
            "local_immersion_mean_A": l_imm["immersion_mean_A"],
            "local_immersion_min_A": l_imm["immersion_min_A"],
            "local_box_center_A_target_frame": l_imm["box_center_A_target_frame"],
        })
        immersed = (l_imm["frac_contact"] >= MIN_FRAC and l_imm["immersion_mean_A"] < MAX_MEAN_A)
        out["local_verdict"] = (
            "local_center_ok" if (l_rmsd <= LOCAL_RMSD_MAX_A and immersed)
            else "local_weak_immersion" if l_rmsd <= LOCAL_RMSD_MAX_A
            else "local_bad_fit")
    else:
        out["local_verdict"] = "too_few_pocket_ca"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default="SUB04")
    ap.add_argument("--pocket-r", type=float, default=12.0)
    ap.add_argument("--json", default="")
    args = ap.parse_args()
    DL_DIR.mkdir(exist_ok=True)

    print(f"pocket_R={args.pocket_r} A | gate local: RMSD<={LOCAL_RMSD_MAX_A} & "
          f"frac>={MIN_FRAC} & imersao<{MAX_MEAN_A} (imersao = backstop honesto)\n")
    hdr = (f"{'ALVO':8}{'SITIO':22}{'ref':6}{'gRMSD':>7}{'gImm':>6}{'gFrc':>6}"
           f"{'nLoc':>5}{'lRMSD':>7}{'lImm':>6}{'lFrc':>6}  local_verdict")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for tid, sname, site, cif in _iter_derived(args.campaign):
        rec = {"target": tid, "site": sname, "ref": site.get("reference_pdb")}
        if not cif or not Path(cif).exists():
            rec["error"] = "no_cif"; rows.append(rec)
            print(f"{tid:8}{sname:22}{str(rec['ref']):6}  ERRO: .cif ausente")
            continue
        try:
            d = diagnose_site(cif, site["reference_pdb"],
                              site.get("reference_ligand_resname", ""), args.pocket_r)
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {str(exc)[:50]}"; rows.append(rec)
            print(f"{tid:8}{sname:22}{str(rec['ref']):6}  ERRO: {rec['error']}")
            continue
        rec.update(d)
        rows.append(rec)
        lr = d.get("local_rmsd_A", "")
        li = d.get("local_immersion_mean_A", "")
        lf = d.get("local_frac_contact", "")
        print(f"{tid:8}{sname:22}{str(d['reference_pdb']):6}"
              f"{d['global_rmsd_A']:>7}{d['global_immersion_mean_A']:>6}"
              f"{d['global_frac_contact']:>6}{d['n_ca_local']:>5}"
              f"{str(lr):>7}{str(li):>6}{str(lf):>6}  {d.get('local_verdict','')}")

    # sumario
    resc = [r for r in rows if r.get("local_verdict") == "local_center_ok"
            and r.get("global_rmsd_A", 0) > 3.0]
    print(f"\nCandidatos a RESGATE por superposicao local (global>3A mas local_center_ok): "
          f"{len(resc)}")
    for r in resc:
        print(f"  - {r['target']} x {r['site']}: global {r['global_rmsd_A']}A -> "
              f"local {r['local_rmsd_A']}A ({r['n_ca_local']} CAs do bolso), "
              f"imersao {r['local_immersion_mean_A']}A frac {r['local_frac_contact']}")
    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"\nJSON -> {args.json}")
    print("\nDIAGNOSTICO read-only: nada escrito no registry.")


if __name__ == "__main__":
    main()
