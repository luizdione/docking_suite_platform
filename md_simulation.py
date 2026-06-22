#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
md_simulation.py - Pipeline OpenMM para TryR + FAD + 2-SCN-NQ (naftoquinona)
============================================================================
Sucessor de md_module.py -> md_clean.py -> md_simulation.py (renomeado em 2026-05-25).
Pipeline end-to-end: setup -> minimizacao -> warm-up multi-fase (0/0.5/1) ->
NVT 50 ps -> NPT 500 ps com release de restraints (k=4->0) -> producao Andersen+Verlet.
Inclui analise pos-rodada (RMSD/RMSF + PBC unwrap do ligante + tempo explicito).

Baseado em:
  * tdudgeon/simple-simulate-complex (MIT)
  * openmm/openmm_workshops section_1
  * openmmforcefields GAFF/SMIRNOFFTemplateGenerator
  * OpenMMDL (JCIM 2024, doi:10.1021/acs.jcim.4c02158)

Uso recomendado: NAO chamar diretamente. Usar o wrapper `run_md_safe.ps1`
(PowerShell como Admin) que aplica o gate `validate_ff.py` + power-cap da GPU
(`nvidia-smi -pl 120`) + serializacao com RNA-seq antes da producao.
"""
from __future__ import annotations
import argparse, logging, math, os, sys, time
from pathlib import Path

try:
    from openmm import (app, unit, LangevinMiddleIntegrator, VerletIntegrator, AndersenThermostat, MonteCarloBarostat, Platform)
    from openmm import CustomExternalForce
except ImportError:
    print("[ERRO] OpenMM nao encontrado. Ative env: conda activate datascience")
    sys.exit(2)

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "docking_project" / "figures" / "md_clean"
WORK.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(WORK / "md_clean_inner.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("md_clean")


# Mapeamento de alvos validados (PDB + n_chains)
TARGET_MAP = {
    "TcTR": {"pdb": "1BZL_clean.pdb", "n_fads": 2, "organism": "T. cruzi"},
    "LiTR": {"pdb": "2JK6_clean.pdb", "n_fads": 2, "organism": "L. infantum"},
    "hGR":  {"pdb": "3DJJ_clean.pdb", "n_fads": 1, "organism": "H. sapiens"},  # 3GRS ausente; 3DJJ tem FAD + ditiol Cys58/Cys63
}


def build_canonical_fad_topology(canonical_template_pdb: Path,
                                 xtal_heavy_pdb: Path,
                                 log_tag: str = "FAD"):
    """Constroi (topology, positions) com nomes canonicos que casam com o XML.

    Estrategia (fix 2026-05-13):
      - `canonical_template_pdb` traz a topologia + nomes corretos (86 atomos
        = 53 heavy + 33 H, identicos ao cofactor_FAD_TcTR.xml).
      - `xtal_heavy_pdb` traz so 53 atomos heavy CRISTALOGRAFICOS, mesmos nomes,
        mesma ordem do template.
      - Substituimos as 53 posicoes heavy pelas do xtal, e propagamos as 33 H
        via transformacao rigida Kabsch (template_heavy -> xtal_heavy).

    Isso evita o problema dos cofactor_FAD_*_xtal_protonated_fixed.pdb, que
    foram renomeados (P1/O1/C1...) e nao casavam com o template do XML
    (PA/O1A/PA...), causando o crash addHydrogens -> No template found.
    """
    import numpy as np
    log.info(f"  [canonical] template={canonical_template_pdb.name} | "
             f"xtal={xtal_heavy_pdb.name}")
    tmpl = app.PDBFile(str(canonical_template_pdb))
    xtal = app.PDBFile(str(xtal_heavy_pdb))

    tmpl_atoms = list(tmpl.topology.atoms())
    xtal_atoms = list(xtal.topology.atoms())
    tmpl_names = [a.name for a in tmpl_atoms]
    xtal_names = [a.name for a in xtal_atoms]
    n_t, n_x = len(tmpl_atoms), len(xtal_atoms)
    n_heavy_tmpl = sum(1 for a in tmpl_atoms if a.element is not None
                       and a.element.symbol != "H")
    n_h_tmpl = n_t - n_heavy_tmpl

    if n_t != 86:
        raise RuntimeError(
            f"{log_tag}: template canonico esperado com 86 atomos "
            f"(53 heavy + 33 H), encontrado {n_t}. "
            f"Verifique {canonical_template_pdb}.")
    if n_x != 53:
        raise RuntimeError(
            f"{log_tag}: xtal heavy esperado com 53 atomos, encontrado {n_x}. "
            f"Verifique {xtal_heavy_pdb}.")

    # Mapa nome -> indice nos dois arquivos
    name_to_tidx = {name: i for i, name in enumerate(tmpl_names)}
    # Os heavy do xtal devem todos existir no template
    missing = [n for n in xtal_names if n not in name_to_tidx]
    if missing:
        raise RuntimeError(
            f"{log_tag}: {len(missing)} atom names do xtal nao existem no template: "
            f"{missing[:8]}{'...' if len(missing)>8 else ''}")

    # Posicoes em nanometros (sem unidades para numpy)
    tmpl_pos = np.array([list(p.value_in_unit(unit.nanometer))
                          for p in tmpl.positions], dtype=float)  # (86,3)
    xtal_pos = np.array([list(p.value_in_unit(unit.nanometer))
                          for p in xtal.positions], dtype=float)  # (53,3)

    # Reordena xtal_pos para casar com os indices do template
    tmpl_heavy_idx = [name_to_tidx[name] for name in xtal_names]  # 53 idx no template
    P = tmpl_pos[tmpl_heavy_idx]   # template heavy positions, ordem do xtal
    Q = xtal_pos                    # xtal heavy positions

    # Kabsch: encontra rotacao R + translacao t que aproxima P @ R + t = Q
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = (Vt.T @ D @ U.T)
    t = cQ - cP @ R
    # RMSD diagnostico (qualidade do match)
    P_aligned = P @ R + t
    rmsd = float(np.sqrt(np.mean(np.sum((P_aligned - Q) ** 2, axis=1))))
    log.info(f"  [canonical] {log_tag}: Kabsch RMSD heavy = {rmsd*10:.3f} A "
             f"(n_heavy={n_x}, n_H_propagated={n_h_tmpl})")

    # ─── Placement dos atomos no frame do xtal ──────────────────────────────
    # FIX 2026-05-13 v2: Kabsch global aplica rotacao rigida unica, MAS o FAD
    # nao e rigido: o template tem uma conformacao genérica e os xtal tem
    # conformacoes dobradas distintas (RMSD 6-9 A entre eles). Propagar Hs
    # via rotacao global coloca-os a 5-10 A do parent heavy correto -> NaN.
    #
    # Solucao: para cada H, ancorar localmente no parent heavy (atomo ao
    # qual esta ligado via Bond na topologia). Vetor H->parent eh rotacionado
    # via R global (que cuida da orientacao macro), mas posicionado no
    # parent_xtal real. Resultado: comprimento de ligacao C-H/N-H/O-H
    # preservado exato, geometria local de grupos rigidos (CH3, CH2,
    # NH, OH) preservada, conformacao macro do FAD = xtal.
    new_pos = np.empty_like(tmpl_pos)
    # 1) Heavy atoms: usar xtal posicao exata
    for k, tidx in enumerate(tmpl_heavy_idx):
        new_pos[tidx] = xtal_pos[k]
    # 2) H atoms: descobrir parent via topology.bonds() e ancorar localmente
    # Indice template -> bonded heavy template indice
    h_indices = [i for i, a in enumerate(tmpl_atoms)
                 if a.element is not None and a.element.symbol == "H"]
    bonds_tmpl = list(tmpl.topology.bonds())
    parent_of = {}  # H tmpl idx -> heavy tmpl idx
    for b in bonds_tmpl:
        i1, i2 = b[0].index, b[1].index
        e1 = tmpl_atoms[i1].element
        e2 = tmpl_atoms[i2].element
        if e1 is not None and e1.symbol == "H" and e2 is not None and e2.symbol != "H":
            parent_of[i1] = i2
        elif e2 is not None and e2.symbol == "H" and e1 is not None and e1.symbol != "H":
            parent_of[i2] = i1
    n_h_ok = 0
    n_h_fallback = 0
    for h_idx in h_indices:
        if h_idx in parent_of:
            p_idx = parent_of[h_idx]              # template heavy idx
            v_tmpl = tmpl_pos[h_idx] - tmpl_pos[p_idx]  # vetor H->parent template
            v_xtal = v_tmpl @ R                   # rotaciona pro frame do xtal
            new_pos[h_idx] = new_pos[p_idx] + v_xtal
            n_h_ok += 1
        else:
            # Fallback: H sem parent claro (improvavel); usar Kabsch global
            new_pos[h_idx] = tmpl_pos[h_idx] @ R + t
            n_h_fallback += 1
    log.info(f"  [canonical] {log_tag}: H placement: {n_h_ok} via parent-anchor, "
             f"{n_h_fallback} via global Kabsch fallback")

    # Reconstroi positions como Quantity nm e devolve a topology do template.
    # Vec3 esta no top-level do openmm 8.x; fallback para openmm.vec3 (7.x).
    try:
        from openmm import Vec3
    except ImportError:
        from openmm.vec3 import Vec3
    positions = [Vec3(*row) * unit.nanometer for row in new_pos]
    return tmpl.topology, positions


def build_system_fast(complex_pdb: Path, fad_xml: Path, lig_xml: Path):
    log.info(f"[FAST] Carregando {complex_pdb.name}")
    pdb = app.PDBFile(str(complex_pdb))
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml",
                        str(fad_xml), str(lig_xml))
    # Inventario do sistema (confirma FAD + LIG + proteina + solvente)
    inv = {}
    for r in pdb.topology.residues():
        inv[r.name] = inv.get(r.name, 0) + 1
    n_atoms = sum(1 for _ in pdb.topology.atoms())
    log.info(f"Sistema: {n_atoms:,} atomos, {sum(inv.values()):,} residuos")
    log.info(f"  proteina (20 AAs): {sum(inv.get(aa,0) for aa in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'])} residuos")
    log.info(f"  cofator FAD:       {inv.get('FAD',0)} molecula(s)")
    log.info(f"  ligante LIG:       {inv.get('LIG',0)} molecula(s)")
    log.info(f"  agua HOH:          {inv.get('HOH',0):,} moleculas")
    log.info(f"  ions Na/Cl:        {inv.get('NA',0)+inv.get('CL',0)}")
    if inv.get('FAD',0) == 0:
        log.warning("  >>> FAD AUSENTE no PDB! Verifique o complex_<TARGET>_initial.pdb")
    log.info("Criando System (PME, 0.9 nm cutoff, HBonds, HMR=4.0 amu)...")
    # [PERF 2026-05-30] HMR (hydrogenMass=4.0 amu) + dt=4 fs na producao -> 2x speedup.
    # Cutoff 0.9 nm (canonico p/ PME) -> ~10% mais rapido sem perda de qualidade.
    # H=4.0 amu e' o valor padrao OpenMM para dt=4fs (H=1.5 amu so e' seguro ate dt=3fs).
    system = ff.createSystem(
        pdb.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        hydrogenMass=4.0 * unit.amu,
        removeCMMotion=False,   # CMMotionRemover pode causar NaN em baixa T
    )
    return pdb, ff, system


def build_system_fresh(receptor_pdb: Path, fad_pdbs: list, lig_pdb: Path,
                       fad_xml: Path, lig_xml: Path):
    """Reconstroi o complexo do ZERO garantindo presenca do FAD.
    1) PDBFixer no receptor (remove ALL HETATM incluindo FAD)
    2) Modeller + cada FAD do fad_pdbs
    3) Modeller + ligante
    4) addHydrogens + addSolvent
    Retorna (Modeller (topology+positions), ForceField, System).
    """
    from pdbfixer import PDBFixer
    log.info(f"[FRESH] Receptor: {receptor_pdb.name}")
    fixer = PDBFixer(filename=str(receptor_pdb))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for k in list(fixer.missingResidues.keys()):
        chain = chains[k[0]]
        if k[1] == 0 or k[1] == len(list(chain.residues())):
            del fixer.missingResidues[k]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # NOTE: removeHeterogens removeria FAD junto; em vez disso, deletar
    # so HOH/GCG/etc manualmente, preservando FAD.
    # FAD cristalografico vem sem H — removemos para usar os pre-protonados
    # (cofactor_FAD_*_xtal_protonated_fixed.pdb, 86 atomos cada).
    HET_TO_DEL = ("HOH", "WAT", "GCG", "NAG", "SO4", "PO4", "EDO", "GOL", "FAD")
    to_del = []
    for res in fixer.topology.residues():
        if res.name in HET_TO_DEL:
            to_del.append(res)
    if to_del:
        log.info(f"  removendo {len(to_del)} HETATM nao-FAD: "
                 f"{sorted(set(r.name for r in to_del))}")
        # PDBFixer delete via Modeller proxy
        from openmm.app import Modeller as _M
        _m = _M(fixer.topology, fixer.positions)
        _m.delete(to_del)
        fixer.topology = _m.topology
        fixer.positions = _m.positions
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # Inventario rapido apos limpeza
    inv_pre = {}
    for r in fixer.topology.residues():
        inv_pre[r.name] = inv_pre.get(r.name, 0) + 1
    log.info(f"  apos PDBFixer: {fixer.topology.getNumResidues()} residuos "
             f"(FAD={inv_pre.get('FAD',0)}, GCG={inv_pre.get('GCG',0)}, "
             f"HOH={inv_pre.get('HOH',0)})")
    log.info(f"  receptor limpo: {fixer.topology.getNumResidues()} residuos, "
             f"{fixer.topology.getNumAtoms()} atomos")

    # ForceField (precisa antes do addHydrogens/addSolvent)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml",
                        str(fad_xml), str(lig_xml))

    # ─────────────────────────────────────────────────────────────────────────
    # ORDEM CRITICA (fix 2026-05-13):
    #
    # Modeller.addHydrogens(ff, ...) chama internamente ff.createSystem() na
    # topologia inteira para validar templates de TODOS os residuos.
    # Se rodarmos addHydrogens DEPOIS de adicionar o FAD pre-protonado,
    # OpenMM tenta casar o FAD do PDB com o template do cofactor_FAD_*.xml
    # — mas os nomes de atomo divergem (PDB: P1/O1/C1...; XML: PA/O1A/PA...)
    # e a chamada explode com:
    #   ValueError: No template found for residue 970 (FAD)
    #
    # Solucao: rodar addHydrogens APENAS na proteina (antes de adicionar
    # FAD/LIG), depois injetar os cofatores ja protonados, e so entao
    # solvatar. Como FAD/LIG so aparecem APOS addHydrogens, eles nunca
    # passam pelo matcher de templates do addHydrogens.
    # O createSystem final (com PME) faz o match dos templates do FAD/LIG
    # a partir dos XMLs correspondentes, na ordem ja existente desses
    # residuos na topologia — sem comparar nomes contra o hydrogens.xml.
    # ─────────────────────────────────────────────────────────────────────────

    # 1) Hidrogenios SO na proteina (sem FAD/LIG ainda na topologia)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    log.info("  addHydrogens(pH=7.4) [protein-only]...")
    modeller.addHydrogens(ff, pH=7.4)
    log.info(f"  proteina protonada: {modeller.topology.getNumAtoms()} atomos "
             f"(was 7459 heavy)")

    # 2) Adicionar FADs pre-protonados (53 heavy + 33 H = 86 atomos canonicos;
    #    Heuristica: se o PDB de entrada tem 53 atomos (heavy-only xtal),
    #    canonicalizar via Kabsch align com cofactor_FAD_TcTR.pdb (86 atomos
    #    canonicos). Se vier ja com 86 atomos canonicos, usar direto.
    # Template PDB canonico. IMPORTANTE: template e fad_xml precisam ser um PAR
    # CONSISTENTE (mesmos elementos), senao o matching quebra. Por isso o v3 PDB
    # so e usado quando o v3 XML tambem existe (mesma logica na resolucao de fad_xml).
    if fad_pdbs:
        _tdir = fad_pdbs[0].parent
        _v3xml = _tdir / "cofactor_FAD_TcTR_v3.xml"
        _v3pdb = _tdir / "cofactor_FAD_TcTR_v3.pdb"
        canonical_template = (_v3pdb if (_v3xml.exists() and _v3pdb.exists())
                              else _tdir / "cofactor_FAD_TcTR.pdb")
    else:
        canonical_template = None
    n_added = 0
    for i, fad_pdb in enumerate(fad_pdbs):
        if not fad_pdb.exists():
            log.warning(f"  FAD #{i+1} ausente: {fad_pdb}; pulando")
            continue
        n_atoms_in = sum(1 for l in open(fad_pdb, "r")
                         if l.startswith(("ATOM  ", "HETATM")))
        if n_atoms_in == 53 and canonical_template and canonical_template.exists():
            log.info(f"  FAD #{i+1}: xtal heavy-only ({n_atoms_in} atomos)"
                     f"; canonicalizando via Kabsch...")
            fad_topo, fad_pos = build_canonical_fad_topology(
                canonical_template, fad_pdb, log_tag=f"FAD#{i+1}")
            modeller.add(fad_topo, fad_pos)
            log.info(f"  FAD #{i+1} adicionado (canonico): "
                     f"{fad_topo.getNumAtoms()} atomos via template {canonical_template.name}")
        else:
            fad = app.PDBFile(str(fad_pdb))
            modeller.add(fad.topology, fad.positions)
            log.info(f"  FAD #{i+1} adicionado (direto): "
                     f"{fad.topology.getNumAtoms()} atomos ({fad_pdb.name})")
        n_added += 1
    if n_added == 0:
        raise RuntimeError(
            f"Nenhum FAD encontrado para {receptor_pdb.name}. "
            f"Para TcTR, espere cofactor_FAD_*_xtal.pdb (53 atomos heavy canonicos) "
            f"+ cofactor_FAD_TcTR.pdb (template canonico 86 atomos) em smoke_test_fad/. "
            f"Para LiTR/hGR, preparar com PDBFixer + addHydrogens primeiro."
        )

    # 3) Adicionar ligante (ja com Hs do lig_TcTR.pdb)
    if lig_pdb.exists():
        lig = app.PDBFile(str(lig_pdb))
        modeller.add(lig.topology, lig.positions)
        log.info(f"  ligante adicionado: {lig.topology.getNumAtoms()} atomos")
    else:
        log.warning(f"  ligante ausente: {lig_pdb}")

    # Solvatacao
    log.info("  addSolvent (TIP3P, padding 1.0 nm, 0.15 M NaCl)...")
    modeller.addSolvent(ff, model="tip3p",
                        padding=1.0 * unit.nanometer,
                        ionicStrength=0.15 * unit.molar,
                        neutralize=True)
    log.info(f"  sistema final: {modeller.topology.getNumAtoms():,} atomos")

    # Inventario
    inv = {}
    for r in modeller.topology.residues():
        inv[r.name] = inv.get(r.name, 0) + 1
    log.info(f"  cofator FAD: {inv.get('FAD',0)}  |  ligante LIG: {inv.get('LIG',0)}  "
             f"|  HOH: {inv.get('HOH',0):,}  |  Na/Cl: {inv.get('NA',0)+inv.get('CL',0)}")
    if inv.get('FAD', 0) == 0:
        raise RuntimeError("FAD nao foi incorporado ao Modeller! Abortando.")

    # [PERF 2026-05-30] cutoff 0.9 nm + HMR=4.0 amu (par com dt=4 fs na producao)
    system = ff.createSystem(
        modeller.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
        hydrogenMass=4.0 * unit.amu,
        removeCMMotion=False,
    )

    # Salvar o PDB reconstruido para inspecao
    out_pdb = WORK / "complex_fresh.pdb"
    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    log.info(f"  PDB reconstruido salvo: {out_pdb.name}")

    # Retorna objeto similar a PDBFile (precisa .topology e .positions)
    class _Wrap:
        def __init__(self, top, pos):
            self.topology = top
            self.positions = pos
    return _Wrap(modeller.topology, modeller.positions), ff, system


def pick_platform(prefer: str = "CUDA", strict: bool = True):
    order = [prefer] if prefer != "AUTO" else ["CUDA", "OpenCL", "CPU"]
    if not strict and prefer not in ("AUTO", None):
        order = [prefer, "OpenCL", "CPU"]

    log.info(f"CONDA_PREFIX={os.environ.get('CONDA_PREFIX','(none)')}")
    log.info(f"OPENMM_PLUGIN_DIR={os.environ.get('OPENMM_PLUGIN_DIR','(none)')}")
    # Probe DLLs Windows
    if sys.platform == "win32":
        cp = os.environ.get("CONDA_PREFIX", "")
        if cp:
            libbin = Path(cp) / "Library" / "bin"
            cudart = list(libbin.glob("cudart64_*.dll")) if libbin.exists() else []
            log.info(f"cudart64_*.dll em {libbin}: "
                     f"{[p.name for p in cudart] if cudart else 'NENHUM'}")
    try:
        from openmm import version as _ommv
        log.info(f"OpenMM version: {_ommv.short_version}")
        try:
            from openmm import Platform as _Plat
            loaded = _Plat.getPluginLoadFailures()
            if loaded:
                log.info(f"Plugin load failures: {list(loaded)}")
            else:
                log.info("Plugins carregados sem falhas.")
        except Exception:
            pass
    except Exception as e:
        log.warning(f"Diagnostico plugins skip: {e}")

    available = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    log.info(f"Plataformas OpenMM disponiveis: {available}")
    if "CUDA" not in available:
        log.warning("CUDA ausente. Conferir:")
        log.warning("  1) cudart64_*.dll no PATH (<env>\\Library\\bin)")
        log.warning("  2) pacote 'openmm-cuda-*' instalado: conda list openmm")
        log.warning("  3) driver NVIDIA suporta a CUDA-toolkit-version do openmm")

    for name in order:
        if name not in available:
            log.warning(f"Plataforma {name} NAO disponivel.")
            continue
        try:
            p = Platform.getPlatformByName(name)
            if name == "CUDA":
                # NOTA: 'single' resolve NaN em sistemas grandes (177k atomos)
                # onde 'mixed' tem bug em algumas versoes/drivers.
                props = {"Precision": "single", "DeterministicForces": "false"}
            elif name == "OpenCL":
                props = {"Precision": "mixed"}
            else:
                props = {}
            log.info(f"Plataforma SELECIONADA: {name}  (props={props})")
            return p, props
        except Exception as e:
            log.warning(f"Falha em {name}: {e}")
    if strict:
        raise RuntimeError(
            f"Plataforma '{prefer}' indisponivel e strict=True. "
            f"Disponiveis: {available}. Use --allow-fallback ou plataforma AUTO."
        )
    return Platform.getPlatformByName("CPU"), {}


def add_position_restraints(system, topology, positions, k=10.0,
                            exclude_resnames=("HOH", "WAT", "NA", "CL", "K",
                                              "FAD", "LIG")):
    # PBC-safe: periodicdistance() respeita a caixa periodica. A forma antiga
    # 0.5*k*((x-x0)^2+...) usa coordenadas absolutas e pode estourar (NaN) no
    # warm-up quando um atomo restringido fica perto da borda da caixa e "envolve"
    # pela PBC -> (x-x0) ~ tamanho da caixa -> forca gigante. (fix 2026-05-25)
    #
    # v2 (2026-05-25): exclui FAD/LIG dos restraints. Restraints em cofatores e
    # ligantes SMIRNOFF impede relaxacao das forcas intra-moleculares durante
    # o warm-up (LG13/LG14 ~1000 kJ/mol/nm) -> NaN imediato. Protocolo padrao
    # de MD: restringir apenas proteina (backbone+sidechain heavy), cofatores
    # e ligantes livres para acomodar no bolso durante equilibracao.
    force = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", k * 4.184 * 100.0)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    exclude_set = set(r.upper() for r in exclude_resnames)
    n = 0
    n_excl = 0
    for atom in topology.atoms():
        if atom.element is None: continue
        if atom.element.symbol == "H": continue
        if atom.residue.name.upper() in exclude_set:
            n_excl += 1
            continue
        pos = positions[atom.index].value_in_unit(unit.nanometer)
        force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
        n += 1
    system.addForce(force)
    log.info(f"Restraints: {n} atomos pesados (proteina) @ k={k} kcal/mol/A^2 "
             f"(excluidos: {n_excl} atomos de {exclude_set - {'HOH','WAT','NA','CL','K'}})")
    return force


def forensic_residue_check(simulation, topology, resnames=("FAD", "LIG")):
    """Diagnostico forense dos residuos custom (FAD/LIG): comprimentos de
    ligacao internos anomalos, contatos proximos (<1.0 A) e forcas por atomo.
    Roda apos a minimizacao para flagrar geometria ruim ANTES da dinamica."""
    import numpy as _np
    state = simulation.context.getState(getPositions=True, getForces=True,
                                        enforcePeriodicBox=False)
    pos = _np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    forces = _np.array(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    fmag = _np.linalg.norm(forces, axis=1)
    atoms = list(topology.atoms())

    target_idx = [a.index for a in atoms if a.residue.name in resnames]
    if not target_idx:
        log.warning(f"[forensic] nenhum residuo {resnames} encontrado")
        return
    log.info(f"[forensic] checando {len(target_idx)} atomos em {resnames}")

    # 1) comprimentos de ligacao internos anomalos
    n_bad_bond = 0
    for b in topology.bonds():
        i1, i2 = b[0].index, b[1].index
        if i1 in target_idx or i2 in target_idx:
            d = float(_np.linalg.norm(pos[i1] - pos[i2])) * 10.0  # A
            if d < 0.7 or d > 2.0:
                n_bad_bond += 1
                a1, a2 = atoms[i1], atoms[i2]
                log.warning(f"[forensic] BOND ANOMALA {d:.2f} A: "
                            f"{a1.name}/{a1.residue.name}{a1.residue.id} - "
                            f"{a2.name}/{a2.residue.name}{a2.residue.id}")
    # 2) contatos proximos: atomo do alvo vs qualquer outro nao-ligado
    bonded = set()
    for b in topology.bonds():
        bonded.add((b[0].index, b[1].index))
        bonded.add((b[1].index, b[0].index))
    n_clash = 0
    for ti in target_idx:
        d = _np.linalg.norm(pos - pos[ti], axis=1) * 10.0  # A
        d[ti] = 99.0
        for nj in _np.where(d < 1.0)[0]:
            if (ti, int(nj)) in bonded:
                continue
            n_clash += 1
            if n_clash <= 15:
                a1, a2 = atoms[ti], atoms[int(nj)]
                log.warning(f"[forensic] CLASH {d[nj]:.2f} A: "
                            f"{a1.name}/{a1.residue.name}{a1.residue.id} - "
                            f"{a2.name}/{a2.residue.name}{a2.residue.id}")
    # 3) forcas nos atomos-alvo
    tgt_arr = _np.array(target_idx)
    tgt_f = fmag[tgt_arr]
    worst = tgt_arr[_np.argsort(tgt_f)[::-1][:8]]
    log.info(f"[forensic] forca max em {resnames}: {tgt_f.max():.2e} kJ/mol/nm "
             f"(media {tgt_f.mean():.2e})")
    for idx in worst:
        a = atoms[idx]
        log.info(f"[forensic]   atom#{idx} F={fmag[idx]:.2e}  "
                 f"{a.name}/{a.residue.name}{a.residue.id}")
    log.info(f"[forensic] resumo: {n_bad_bond} bonds anomalas, {n_clash} clashes <1.0 A")
    if n_bad_bond or n_clash:
        log.warning("[forensic] >>> geometria ruim em FAD/LIG detectada "
                    "(provavel causa do NaN no warm-up)")



def diagnose_system(simulation, topology, tag="post-min"):
    """Imprime top-10 forcas + box vectors + atomos NaN. Use apos minimizacao."""
    import numpy as _np
    state = simulation.context.getState(getPositions=True, getForces=True,
                                         enforcePeriodicBox=False)
    pos = _np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    forces = _np.array(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    fmag = _np.linalg.norm(forces, axis=1)

    # Box vectors
    box = state.getPeriodicBoxVectors()
    log.info(f"[diag/{tag}] box vectors (nm): a={box[0]} b={box[1]} c={box[2]}")

    # NaN scan
    nan_pos = _np.isnan(pos).any(axis=1)
    nan_force = _np.isnan(forces).any(axis=1)
    log.info(f"[diag/{tag}] atoms with NaN positions: {nan_pos.sum()}")
    log.info(f"[diag/{tag}] atoms with NaN forces:    {nan_force.sum()}")
    if nan_pos.sum() or nan_force.sum():
        bad = _np.where(nan_pos | nan_force)[0][:10]
        atoms = list(topology.atoms())
        for idx in bad:
            a = atoms[idx]
            log.error(f"  NaN at atom #{idx}: {a.name} of {a.residue.name}"
                      f" {a.residue.id} chain {a.residue.chain.id}")

    # Top-10 forcas
    fmag_finite = _np.where(_np.isnan(fmag), 0.0, fmag)
    top_idx = _np.argsort(fmag_finite)[::-1][:10]
    atoms = list(topology.atoms())
    log.info(f"[diag/{tag}] TOP-10 forcas (kJ/mol/nm):")
    for i, idx in enumerate(top_idx):
        a = atoms[idx]
        log.info(f"  #{i+1} atom#{idx} F={fmag_finite[idx]:.2e}"
                 f"  {a.name}/{a.residue.name}{a.residue.id}"
                 f"  chain={a.residue.chain.id}")
    # Energia potencial total
    state_e = simulation.context.getState(getEnergy=True)
    pe = state_e.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    log.info(f"[diag/{tag}] PotentialEnergy = {pe:.2e} kJ/mol")
    return fmag_finite.max(), nan_pos.sum() + nan_force.sum()

def _resume_production(pdb, system, ns_prod: float, tag: str,
                       platform, props, chk_path: Path, mode: str = None):
    """Retoma a PRODUCAO de um checkpoint .chk apos queda de energia, sem refazer
    minimizacao/warm-up/equilibracao. Reconstroi as forcas EXATAMENTE como estavam
    no fim da equilibracao (restraint[k] + Andersen + Barostat) para o checkpoint
    casar, carrega o .chk, e roda apenas os steps restantes (append nos reporters).
    Recupera os ns ja simulados antes da queda."""
    # TRAVA DE MODO: nao deixa retomar um checkpoint de OUTRA hipotese (ex.: pedir
    # 'subversive' mas o .chk ser do run inibidor) -> contaminaria a trajetoria.
    mode_file = WORK / f"md_{tag}.mode"
    if mode and mode_file.exists():
        prev = mode_file.read_text().strip()
        if prev and prev != mode:
            raise RuntimeError(
                f"checkpoint do tag '{tag}' e do modo '{prev}', mas voce pediu '{mode}'. "
                f"Use outro --tag (ex.: TcTR_{mode}) para nao sobrescrever/contaminar o run '{prev}'.")
    use_restraint = os.environ.get("MD_CLEAN_NO_RESTRAINT", "0") != "1"
    if use_restraint:
        add_position_restraints(system, pdb.topology, pdb.positions, k=5.0)
    # [PERF 2026-05-30] LangevinMiddleIntegrator (substitui Verlet+Andersen) + dt=4 fs
    # + Barostat a cada 100 steps (era 25). Sem AndersenThermostat -> menos forcas extras.
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 100))
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                          4.0 * unit.femtosecond)
    simulation = app.Simulation(pdb.topology, system, integrator, platform, props)
    simulation.context.setPositions(pdb.positions)  # placeholder (sobrescrito pelo chk)
    simulation.loadCheckpoint(str(chk_path))         # restaura pos/vel/box/params/step

    done = simulation.currentStep
    # dt=4 fs -> 250.000 steps/ns (era 500.000 com dt=2 fs)
    prod_steps = int(ns_prod * 1000 * unit.picosecond / (4 * unit.femtosecond))
    remaining = prod_steps - done
    log.info(f"[resume] checkpoint em step {done:,}/{prod_steps:,} "
             f"({100.0*done/max(prod_steps,1):.1f}% | {done*4/1e6:.2f} ns ja feitos).")
    dcd_path = WORK / f"md_{tag}.dcd"
    log_path = WORK / f"md_{tag}.csv"
    if remaining <= 0:
        log.info("[resume] producao ja estava completa; nada a fazer.")
        return dcd_path, log_path

    # Alinha o DCD ao checkpoint (descarta frames pos-checkpoint que ficaram orfaos
    # da queda) para nao duplicar/desordenar frames ao dar append.
    try:
        import mdtraj as _md
        topo = WORK / f"minimized_{tag}.pdb"
        if topo.exists() and dcd_path.exists():
            keep = done // 2500  # DCDReporter grava a cada 2500 steps (= 10 ps)
            tr = _md.load(str(dcd_path), top=str(topo))
            if tr.n_frames > keep:
                tr[:keep].save_dcd(str(dcd_path))
                log.info(f"[resume] DCD truncado p/ {keep} frames (alinhado ao checkpoint).")
    except Exception as e:
        log.warning(f"[resume] truncagem do DCD pulada ({e}); pode haver leve sobreposicao de frames.")

    # [PERF] DCD a cada 2500 steps = 10 ps (mesma cadencia que antes); CSV a cada 25.000 steps
    simulation.reporters.append(app.DCDReporter(str(dcd_path), 2500, append=True))
    simulation.reporters.append(app.StateDataReporter(
        str(log_path), 25000, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=prod_steps, separator=",", append=True))
    simulation.reporters.append(app.CheckpointReporter(str(chk_path), 25000))
    log.info(f"[resume] continuando producao por {remaining:,} steps...")
    t0 = time.time()
    simulation.step(remaining)
    el = time.time() - t0
    if el > 0:
        log.info(f"  Producao (resume) OK em {el/60:.1f} min")
    return dcd_path, log_path


def run_md(pdb, ff, system, ns_prod: float, tag: str = "TcTR",
           platform_name: str = "CUDA", strict_platform: bool = True,
           resume: bool = True, mode: str = None):
    platform, props = pick_platform(prefer=platform_name, strict=strict_platform)

    # ── RESUME pos-queda: se ha checkpoint, retoma a producao de onde parou ──
    chk_path = WORK / f"md_{tag}.chk"
    if resume and chk_path.exists():
        log.info(f"[resume] checkpoint encontrado ({chk_path.name}) -> tentando retomar.")
        try:
            return _resume_production(pdb, system, ns_prod, tag, platform, props, chk_path, mode)
        except Exception as e:
            log.error(f"[resume] FALHOU ({e}). Para recomecar do zero, rode com --no-resume "
                      f"(isso ignora/sobrescreve o checkpoint). Abortando para nao corromper.")
            raise
    # Restraint na PROTEINA (exclui FAD/LIG — v2 fix NaN warm-up).
    # MD_CLEAN_NO_RESTRAINT=1 desliga tudo (debug only).
    use_restraint = os.environ.get("MD_CLEAN_NO_RESTRAINT", "0") != "1"
    if use_restraint:
        add_position_restraints(system, pdb.topology, pdb.positions, k=5.0)
    else:
        log.warning("[diag] Restraints DESLIGADOS via MD_CLEAN_NO_RESTRAINT=1")
    # [PERF 2026-05-30] LangevinMiddleIntegrator no warmup E producao (substitui
    # Verlet + AndersenThermostat -> ~15% overhead a menos, mais estavel com HMR).
    # NOTA HISTORICA: o Langevin "classico" causava NaN nesse sistema 177k. O
    # LangevinMiddle (split-operator) usa update OBABO e estabiliza nucleo PME
    # melhor. Combinado com HMR=4 amu + dt=4 fs, e' o protocolo canonico OpenMM 8.x.
    integrator = LangevinMiddleIntegrator(100 * unit.kelvin, 1.0 / unit.picosecond,
                                          1.0 * unit.femtosecond)  # dt=1 fs warm-up
    simulation = app.Simulation(pdb.topology, system, integrator, platform, props)
    simulation.context.setPositions(pdb.positions)
    if pdb.topology.getPeriodicBoxVectors() is not None:
        simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

    # ── (1) Minimizacao apertada ─────────────────────────────────────────────
    log.info("[1/5] Minimizacao apertada (max 20000 iter, tol=10 kJ/mol/nm)...")
    t0 = time.time()
    simulation.minimizeEnergy(
        tolerance=10.0 * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=20000,
    )
    log.info(f"  ok ({time.time()-t0:.1f} s)")
    # Salva snapshot pos-minimizacao p/ debug
    try:
        state = simulation.context.getState(getPositions=True)
        with open(WORK / f"minimized_{tag}.pdb", "w") as fh:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), fh)
    except Exception as e:
        log.warning(f"snapshot post-min falhou: {e}")

    # ── DIAGNOSTICO forense pos-minimizacao ──────────────────────────────
    max_f, n_nan = diagnose_system(simulation, simulation.topology, tag="post-min")
    if n_nan > 0:
        raise RuntimeError(f"Sistema ja tem {n_nan} NaN apos minimizacao — XMLs/PDB corrompidos.")
    if max_f > 1e7:
        log.warning(f"[diag] Forca maxima {max_f:.2e} kJ/mol/nm > 1e7 — clash provavel.")
        log.warning("[diag] Tentando minimizacao extra com tol=1.0 kJ/mol/nm...")
        simulation.minimizeEnergy(
            tolerance=1.0 * unit.kilojoule_per_mole / unit.nanometer,
            maxIterations=50000,
        )
        max_f, n_nan = diagnose_system(simulation, simulation.topology, tag="post-min-2")
        if max_f > 1e7 or n_nan > 0:
            raise RuntimeError(f"Forca ainda {max_f:.2e}, NaN={n_nan}. Inspecione minimized_{tag}.pdb")

    # ── (2) Warm-up gradual 50→100→200→300 K (NVT, dt=1 fs) ───────────────────
    log.info("[2/5] Warm-up: diagnostico forense + fase 0 conservadora")
    # ── DIAGNOSTICO FORENSE de FAD/LIG (residuos custom) ──
    forensic_residue_check(simulation, simulation.topology, ("FAD", "LIG"))

    # helper: detecta NaN e aborta apontando o 1o atomo ruim
    def _nan_guard(stage):
        import numpy as _np
        st = simulation.context.getState(getPositions=True)
        p = _np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        mask = _np.isnan(p).any(axis=1)
        nbad = int(mask.sum())
        if nbad:
            atoms = list(simulation.topology.atoms())
            for idx in _np.where(mask)[0][:12]:
                a = atoms[int(idx)]
                log.error(f"  [{stage}] NaN @ atom#{idx} {a.name}/"
                          f"{a.residue.name}{a.residue.id} chain={a.residue.chain.id}")
            raise RuntimeError(f"NaN em '{stage}': {nbad} atomos (ver IDs acima). "
                               f"Se FAD/LIG -> geometria custom ruim; "
                               f"se HOH -> clash de solvatacao.")

    # [PERF 2026-05-30] Warmup ENXUTO (fases 0 e 0.5 removidas; comeca direto em 100 K).
    # Justificativa: sistema vem de minimizacao apertada (tol=10 kJ/mol/nm); sem clashes
    # severos, dt=0.25fs @ 5K era over-conservador (gastava 1 h por run). Langevin Middle
    # acomoda starting kinetic energy bem em 100 K com dt=1 fs.
    log.info("[2/5] Warm-up gradual NVT 100->200->300 K (Langevin, dt=1 fs)...")
    integrator.setStepSize(1.0 * unit.femtosecond)
    for T in (100, 200, 300):
        integrator.setTemperature(T * unit.kelvin)
        simulation.context.setVelocitiesToTemperature(T * unit.kelvin)
        n = int(25 * unit.picosecond / (1 * unit.femtosecond))  # 25 ps cada
        simulation.step(n)
        _nan_guard(f"warmup T={T}K")
        log.info(f"  warm-up T={T} K ({n:,} steps de 1 fs) ok")

    # ── (3) NVT 50 ps consolidacao a 300 K, dt=1 fs ──
    log.info("[3/5] NVT consolidacao 50 ps @ 300 K (dt=1 fs)...")
    simulation.step(int(50 * unit.picosecond / (1 * unit.femtosecond)))
    _nan_guard("consolidacao NVT")
    log.info("  consolidacao NVT ok")

    # ── (4) NPT 500 ps com release gradual do restraint, dt=2 fs ─────────────
    # [PERF 2026-05-30] Barostat a cada 100 steps (era 25) -- aceitancia ainda alta,
    # menos overhead. dt=2 fs no NPT: HMR esta presente mas mantemos 2 fs aqui por
    # margem extra durante o release; 4 fs so e' usado na producao.
    log.info("[4/5] NPT 500 ps @ 1 bar (release gradual restraint, dt=2 fs)...")
    integrator.setStepSize(2.0 * unit.femtosecond)
    integrator.setTemperature(300 * unit.kelvin)
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 100))
    simulation.context.reinitialize(preserveState=True)
    # Release gradual do restraint SO se ele existir (com MD_CLEAN_NO_RESTRAINT=1
    # nao ha parametro 'k' -> faz NPT direto). Evita 'invalid parameter name: k'.
    if "k" in simulation.context.getParameters():
        for k_kcal in (4.0, 2.0, 1.0, 0.5, 0.0):
            simulation.context.setParameter("k", k_kcal * 4.184 * 100.0)
            simulation.step(int(100 * unit.picosecond / (2 * unit.femtosecond)))
            log.info(f"  release k={k_kcal} kcal/mol/A^2 ok")
    else:
        log.info("  (sem restraint) NPT 500 ps direto @ 1 bar...")
        simulation.step(int(500 * unit.picosecond / (2 * unit.femtosecond)))
        log.info("  NPT 500 ps ok")

    # ── (5) Producao dt=4 fs com HMR + Langevin Middle ───────────────────────
    # [PERF 2026-05-30] dt=4 fs (era 2 fs) + HMR=4 amu -> ~2x speedup.
    # DCDReporter a cada 2500 steps = 10 ps/frame (mesma cadencia que antes);
    # StateDataReporter a cada 25000 steps = 100 ps/linha CSV (era 5 ps -> 20x menos
    # I/O); CheckpointReporter a cada 25000 steps = 100 ps (era 100 ps tambem).
    integrator.setStepSize(4.0 * unit.femtosecond)
    integrator.setTemperature(300 * unit.kelvin)
    log.info("[5/5] Producao iniciando (dt=4 fs, Langevin Middle, HMR=4 amu)...")
    prod_steps = int(ns_prod * 1000 * unit.picosecond / (4 * unit.femtosecond))
    dcd_path = WORK / f"md_{tag}.dcd"
    log_path = WORK / f"md_{tag}.csv"
    chk_path = WORK / f"md_{tag}.chk"
    simulation.reporters.append(app.DCDReporter(str(dcd_path), 2500))
    simulation.reporters.append(app.StateDataReporter(
        str(log_path), 25000, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=prod_steps, separator=","))
    simulation.reporters.append(app.CheckpointReporter(str(chk_path), 25000))

    # Sidecar de modo: registra a hipotese deste run p/ a trava de resume
    if mode:
        (WORK / f"md_{tag}.mode").write_text(str(mode))

    log.info(f"  Producao {ns_prod} ns (dt=4 fs, {prod_steps:,} steps)...")
    t0 = time.time()
    simulation.step(prod_steps)
    el = time.time() - t0
    if el > 0:
        log.info(f"  Producao OK em {el/60:.1f} min ({ns_prod*86400/el:.1f} ns/dia)")
    return dcd_path, log_path


def analyze(dcd_path: Path, topology_pdb: Path, tag: str = "TcTR",
            frame_interval_ps: float = 10.0):
    """Analise da trajetoria de MD:
      - RMSD(t) da proteina (backbone CA) com media movel;
      - RMSD do LIGANTE NO BOLSO: medido APOS (i) PBC unwrap minimum-image do
        ligante em relacao ao centro da proteina e (ii) superposicao na proteina,
        SEM re-superpor no proprio ligante -> capta deriva/saida do sitio (e nao
        mudanca de conformacao interna, que era o que a versao antiga media);
      - RMSF por residuo (CA);
      - curva de CONVERGENCIA (media cumulativa do RMSD; plato = equilibrado) e
        estatistica de equilibrio da 2a metade da trajetoria.

    Param `frame_interval_ps`: intervalo real entre frames do DCD (default 10 ps =
    5000 steps * 2 fs, que e o DCDReporter de producao deste pipeline). Usado para
    construir o eixo de tempo PORQUE o `traj.time` do DCD escrito por OpenMM 8.2 +
    plugin CUDA 110 vem com fator de 10 errado neste sistema; nao confiar nele.

    Salva PNG (3 paineis) + resumo JSON.
    """
    try:
        import numpy as np
        import mdtraj as md
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("numpy/mdtraj/matplotlib indisponiveis; pulando analise.")
        return
    log.info("Carregando trajetoria para analise...")
    traj = md.load(str(dcd_path), top=str(topology_pdb))
    if traj.n_frames < 2:
        log.warning(f"  trajetoria com {traj.n_frames} frame(s); analise pulada.")
        return

    ca = traj.topology.select("protein and name CA")
    lig = traj.topology.select("resname LIG")

    # ── PBC unwrap do ligante (minimum-image) ANTES da superposicao ───────────
    # Sem isso, quando o ligante atravessa a borda periodica, traj.xyz salta o
    # comprimento da caixa e o RMSD do ligante explode para 100+ A (artefato puro).
    # Fazemos minimum-image SO no ligante, relativo ao centro de massa da proteina,
    # frame a frame. Funciona aqui porque a caixa e cubica/ortogonal.
    if len(lig) and traj.unitcell_lengths is not None:
        prot_com = traj.xyz[:, ca, :].mean(axis=1)            # (nF, 3) nm
        lig_com  = traj.xyz[:, lig, :].mean(axis=1)           # (nF, 3) nm
        box      = traj.unitcell_lengths                      # (nF, 3) nm
        delta    = lig_com - prot_com                          # (nF, 3)
        shift    = np.round(delta / box) * box                # imagem mais proxima
        traj.xyz[:, lig, :] = traj.xyz[:, lig, :] - shift[:, None, :]
        max_shift_A = float(np.abs(shift).max()) * 10.0
        if max_shift_A > 0:
            log.info(f"  [analise] PBC unwrap ligante: max |shift| = {max_shift_A:.2f} A")

    # Alinhar TODOS os frames na proteina (CA) -> referencial do receptor.
    # traj.xyz passa a estar no referencial da proteina (in-place).
    traj.superpose(traj, 0, atom_indices=ca)

    # Tempo: usa frame_interval_ps explicito (ver docstring sobre traj.time bug).
    time_ns = (np.arange(1, traj.n_frames + 1) * frame_interval_ps) / 1000.0

    # RMSD da proteina (backbone CA)
    rmsd_prot = md.rmsd(traj, traj, 0, atom_indices=ca) * 10.0   # nm -> A

    # RMSD do ligante NO BOLSO (deriva): distancia no referencial ja alinhado na
    # proteina, frame a frame vs frame 0 -- NAO re-superpoe no ligante.
    rmsd_lig_pocket = rmsd_lig_intra = None
    if len(lig):
        d = traj.xyz[:, lig, :] - traj.xyz[0, lig, :]
        rmsd_lig_pocket = np.sqrt((d ** 2).sum(axis=2).mean(axis=1)) * 10.0
        rmsd_lig_intra = md.rmsd(traj, traj, 0, atom_indices=lig) * 10.0  # conformacao

    # Convergencia: media cumulativa e media movel do RMSD da proteina
    n = traj.n_frames
    cummean = np.cumsum(rmsd_prot) / np.arange(1, n + 1)
    win = max(1, n // 20)
    runmean = np.convolve(rmsd_prot, np.ones(win) / win, mode="same")

    # Equilibrio: estatistica da 2a metade
    half = n // 2
    eq_mean = float(rmsd_prot[half:].mean())
    eq_std = float(rmsd_prot[half:].std())

    # RMSF por residuo (CA), na trajetoria ja alinhada
    rmsf = md.rmsf(traj, traj, 0, atom_indices=ca) * 10.0
    res_ids = [traj.topology.atom(int(i)).residue.resSeq for i in ca]

    # ---- figura: 3 paineis ----
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    ax[0].plot(time_ns, rmsd_prot, color="#1f77b4", lw=0.8, alpha=0.55, label="Proteina (CA)")
    ax[0].plot(time_ns, runmean, color="#08306b", lw=1.8, label=f"media movel ({win}f)")
    if rmsd_lig_pocket is not None:
        ax[0].plot(time_ns, rmsd_lig_pocket, color="#d62728", lw=1.0, alpha=0.85,
                   label="Ligante no bolso")
    ax[0].axvspan(float(time_ns[half]), float(time_ns[-1]), color="green", alpha=0.06)
    ax[0].set(xlabel="Tempo (ns)", ylabel="RMSD (A)",
              title=f"RMSD - {tag}  (equil.: {eq_mean:.2f}+/-{eq_std:.2f} A)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    ax[1].plot(res_ids, rmsf, color="#2ca02c")
    ax[1].set(xlabel="Residuo (resSeq)", ylabel="RMSF (A)", title=f"RMSF por residuo - {tag}")
    ax[1].grid(alpha=0.3)

    ax[2].plot(time_ns, cummean, color="#6a3d9a", lw=1.8)
    ax[2].axhline(eq_mean, ls="--", color="gray", lw=1, label=f"media equil. {eq_mean:.2f} A")
    ax[2].set(xlabel="Tempo (ns)", ylabel="RMSD medio cumulativo (A)",
              title=f"Convergencia - {tag}  (plato = equilibrado)")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    fig.tight_layout()
    out_png = WORK / f"md_{tag}_rmsd_rmsf.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # resumo numerico
    import json
    summary = {
        "tag": tag, "n_frames": int(n), "t_total_ns": float(time_ns[-1]),
        "rmsd_prot_final_A": float(rmsd_prot[-1]),
        "rmsd_prot_equil_mean_A": eq_mean, "rmsd_prot_equil_std_A": eq_std,
        "rmsd_lig_pocket_final_A": (None if rmsd_lig_pocket is None else float(rmsd_lig_pocket[-1])),
        "rmsd_lig_pocket_max_A": (None if rmsd_lig_pocket is None else float(rmsd_lig_pocket.max())),
        "rmsd_lig_intra_final_A": (None if rmsd_lig_intra is None else float(rmsd_lig_intra[-1])),
        "rmsf_mean_A": float(rmsf.mean()), "rmsf_max_A": float(rmsf.max()),
    }
    (WORK / f"md_{tag}_rmsd_summary.json").write_text(json.dumps(summary, indent=2))
    log.info(f"  RMSD/RMSF/convergencia salvo em {out_png.name}")
    _ligtxt = "n/a" if rmsd_lig_pocket is None else f"{rmsd_lig_pocket[-1]:.2f} A"
    log.info(f"  [analise] proteina equil={eq_mean:.2f}+/-{eq_std:.2f} A | "
             f"ligante no bolso (final)={_ligtxt} | RMSF max={rmsf.max():.2f} A")


def resolve_topology(args) -> Path:
    """Escolhe o PDB de topologia que CASA com o DCD (mesmo Nde atomos).

    Bug historico: analyze() recebia args.complex (complex_<TARGET>_initial.pdb),
    que NAO existe em rodadas --fresh -> OSError. A topologia correta e a do
    sistema solvatado realmente simulado:
      1. minimized_<tag>.pdb  (escrito por run_md apos a minimizacao)   <- ideal
      2. complex_fresh.pdb    (escrito por build_system_fresh)
      3. args.complex         (modo --fast: ja e o sistema pronto)
    """
    candidates = [
        WORK / f"minimized_{args.tag}.pdb",
        WORK / "complex_fresh.pdb",
        args.complex,
    ]
    for c in candidates:
        if c and Path(c).exists():
            return Path(c)
    # ultimo recurso: devolve args.complex (analyze loga o erro de forma clara)
    return Path(args.complex)


def main():
    p = argparse.ArgumentParser(description="md_clean.py - MD enxuto TryR/FAD/LIG")
    p.add_argument("--target", default="TcTR", choices=list(TARGET_MAP.keys()),
                   help="Alvo (TcTR/LiTR/hGR). Default: TcTR")
    p.add_argument("--ns", type=float, default=5.0, help="ns de producao")
    p.add_argument("--tag", default="TcTR", help="Sufixo dos arquivos de saida")
    p.add_argument("--platform", default="CUDA",
                   choices=["CUDA", "OpenCL", "CPU", "AUTO"],
                   help="Plataforma OpenMM")
    p.add_argument("--allow-fallback", action="store_true",
                   help="Permite cair para CPU/OpenCL se plataforma pedida indisponivel")
    p.add_argument("--complex", type=Path, default=None,
                   help="Path para o complex_<TARGET>_initial.pdb (modo --fast)")
    p.add_argument("--fad-xml", type=Path, default=None)
    p.add_argument("--lig-xml", type=Path, default=None)
    p.add_argument("--receptor", type=Path, default=None,
                   help="Path para o receptor_clean.pdb (modo --fresh)")
    p.add_argument("--fad1", type=Path, default=None,
                   help="FAD #1 PDB pre-protonado canonico (modo --fresh)")
    p.add_argument("--fad2", type=Path, default=None,
                   help="FAD #2 PDB pre-protonado canonico (modo --fresh)")
    p.add_argument("--lig-pdb", type=Path, default=None,
                   help="Ligante PDB (modo --fresh)")
    p.add_argument("--analyze-only", action="store_true",
                   help="So roda analyze() sobre dcd existente")
    p.add_argument("--fast", action="store_true",
                   help="Carrega complex_<TARGET>_initial.pdb ja pronto")
    p.add_argument("--fresh", action="store_true",
                   help="Reconstroi sistema do zero (PDBFixer + Modeller + FAD + LIG)")
    p.add_argument("--mode", default=None, choices=["inhibitor", "subversive"],
                   help="Hipotese de ligacao: 'inhibitor' (classico/inibidor de TR, "
                        "pose periferica lig_<TARGET>.pdb) ou 'subversive' (substrato "
                        "subversivo, pose na si-face do FAD lig_<TARGET>_fadsite.pdb). "
                        "Se omitido e o terminal for interativo, sera perguntado.")
    p.add_argument("--no-resume", action="store_true",
                   help="Ignora checkpoint .chk existente e recomeca do zero "
                        "(default: retoma producao de onde parou apos queda).")
    args = p.parse_args()
    if args.fresh:
        args.fast = False  # FRESH ganha precedencia

    # ── MODO de ligacao (inibidor vs substrato subversivo) ───────────────────
    # Define a pose inicial do ligante e o sufixo do tag. Pergunta se nao vier
    # por --mode e o terminal for interativo (caixa de preenchimento no CLI).
    if args.mode is None and not args.analyze_only:
        if sys.stdin and sys.stdin.isatty():
            print("\n" + "=" * 60)
            print(" MODO DE ANALISE DA MD (TcTR + 2-SCN-NQ)")
            print("   [1] classico / INIBIDOR de TR   (pose periferica)")
            print("   [2] substrato SUBVERSIVO        (si-face do FAD)")
            print("=" * 60)
            _r = input(" Escolha [1/2] (default 1): ").strip()
            args.mode = "subversive" if _r == "2" else "inhibitor"
        else:
            args.mode = "inhibitor"
            log.info("[modo] nao-interativo e sem --mode -> default 'inhibitor'.")
    log.info(f"[modo] {args.mode}")
    # separa as saidas das duas hipoteses quando o tag e o default
    if args.mode and args.tag == "TcTR":
        args.tag = "TcTR_%s" % args.mode

    # TRAVA DE MODO ANTECIPADA (antes do build caro de ~45 s): se vamos retomar
    # e o checkpoint desse tag e de OUTRO modo, aborta JA com instrucao clara.
    if not args.analyze_only and not args.no_resume and args.mode:
        _chk = WORK / f"md_{args.tag}.chk"
        _modef = WORK / f"md_{args.tag}.mode"
        if _chk.exists() and _modef.exists():
            _prev = _modef.read_text().strip()
            if _prev and _prev != args.mode:
                log.error("=" * 70)
                log.error(f"ABORTADO: o checkpoint do tag '{args.tag}' e do modo '{_prev}', "
                          f"mas voce pediu '{args.mode}'.")
                log.error(f"  -> Use um TAG distinto (ex.: --tag TcTR_{args.mode}) para nao "
                          f"sobrescrever/contaminar o run '{_prev}'.")
                log.error(f"  -> Ou rode com --no-resume para recomecar este tag do zero.")
                log.error("=" * 70)
                sys.exit(4)

    log.info("=" * 70)
    log.info("md_clean.py - TryR + FAD + 2-SCN-NQ")
    log.info(f"Workdir: {WORK}")
    log.info(f"Args: ns={args.ns}, tag={args.tag}, platform={args.platform}, "
             f"allow_fallback={args.allow_fallback}, fresh={args.fresh}, fast={args.fast}")
    log.info("=" * 70)

    # Resolve paths default por target (TARGET_MAP)
    base = ROOT / "docking_project/figures"
    proteins_dir = ROOT / "docking_project/proteins"
    if args.target not in TARGET_MAP:
        log.warning(f"Target {args.target} nao mapeado; usando default TcTR.")
        args.target = "TcTR"
    target_info = TARGET_MAP[args.target]
    if args.target == "TcTR":
        sub = base / "smoke_test_fad"
    else:
        sub = base / f"smoke_test_{args.target.lower()}"
    tagged = args.target
    if args.complex is None:
        args.complex = sub / f"complex_{tagged}_initial.pdb"
    if args.fad_xml is None:
        args.fad_xml = sub / f"cofactor_FAD_{tagged}.xml"
        if not args.fad_xml.exists():
            args.fad_xml = base / "smoke_test_fad" / "cofactor_FAD_TcTR.xml"
            log.info(f"  reusando FAD XML do TcTR para {args.target}")
    if args.lig_xml is None:
        args.lig_xml = sub / f"lig_{tagged}_smirnoff.xml"
        if not args.lig_xml.exists():
            args.lig_xml = base / "smoke_test_fad" / "lig_TcTR_smirnoff.xml"
            log.info(f"  reusando ligante XML do TcTR para {args.target}")
    log.info(f"Target: {args.target} ({target_info['organism']}, "
             f"{target_info['n_fads']} FAD(s) esperado(s))")
    log.info(f"  receptor: {target_info['pdb']}")
    log.info(f"  FAD XML:  {args.fad_xml.name}")
    log.info(f"  LIG XML:  {args.lig_xml.name}")

    if args.analyze_only:
        _topo = resolve_topology(args)
        log.info(f"  [analyze-only] topologia: {_topo.name}")
        analyze(WORK / f"md_{args.tag}.dcd", _topo, args.tag)
        return

    # Precondicoes de arquivo
    if args.fast and not args.complex.exists():
        log.warning(f"[FAST] {args.complex.name} nao existe - caindo para FRESH automaticamente.")
        args.fast = False
        args.fresh = True

    if args.fresh:
        # FRESH: precisa apenas dos XMLs (PDB do receptor vem do TARGET_MAP)
        for f in (args.fad_xml, args.lig_xml):
            if not f.exists():
                log.error(f"Arquivo ausente: {f}")
                sys.exit(3)
    else:
        for f in (args.complex, args.fad_xml, args.lig_xml):
            if not f.exists():
                log.error(f"Arquivo ausente: {f}")
                sys.exit(3)

    if args.fresh:
        receptor = args.receptor or (proteins_dir / TARGET_MAP[args.target]["pdb"])
        # NOTA (fix 2026-05-13): por default agora usamos os xtal heavy PDBs
        # (53 atomos canonicos) + template TcTR.pdb (86 atomos canonicos), que
        # batem com cofactor_FAD_TcTR.xml. Os antigos
        # cofactor_FAD_*_xtal_protonated_fixed.pdb (88 atomos renomeados)
        # NAO casam com o XML e causam o crash addHydrogens. O helper
        # build_canonical_fad_topology() faz Kabsch align template->xtal.
        fad1 = args.fad1 or (sub / "cofactor_FAD_0_xtal.pdb")
        fad2 = args.fad2 or (sub / "cofactor_FAD_1_xtal.pdb")
        # MODO define a pose inicial do ligante (subversivo = na si-face do FAD)
        if args.lig_pdb:
            lig_pdb = args.lig_pdb
        elif args.mode == "subversive":
            lig_pdb = sub / f"lig_{tagged}_fadsite.pdb"
        else:
            lig_pdb = sub / f"lig_{tagged}.pdb"
        if not lig_pdb.exists():
            log.error(f"Pose inicial do ligante ausente para modo '{args.mode}': {lig_pdb}")
            if args.mode == "subversive":
                log.error("  Gere com build_fad_siface_pose.py (kit fig13_3d_kit).")
            sys.exit(3)
        log.info(f"[FRESH] receptor={receptor.name} fad1={fad1.name} fad2={fad2.name} lig={lig_pdb.name}")
        pdb, ff, system = build_system_fresh(receptor, [fad1, fad2], lig_pdb,
                                              args.fad_xml, args.lig_xml)
    else:
        pdb, ff, system = build_system_fast(args.complex, args.fad_xml, args.lig_xml)
    dcd_path, _ = run_md(pdb, ff, system, args.ns, args.tag,
                         platform_name=args.platform,
                         strict_platform=(not args.allow_fallback),
                         resume=(not args.no_resume), mode=args.mode)
    _topo = resolve_topology(args)
    log.info(f"  [analyze] topologia: {_topo.name}")
    analyze(dcd_path, _topo, args.tag)
    log.info("DONE.")


if __name__ == "__main__":
    main()
