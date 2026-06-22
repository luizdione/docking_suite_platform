#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replica_driver.py - Roda REPLICAS INDEPENDENTES de MD das poses subversivas
(TcTR, LiTR, hGR) sequencialmente, robusto a quedas, por ~3 dias sem supervisao.

Motivacao cientifica
---------------------
Responder a critica de revisor sobre "trajetoria unica (n=1)" nas MDs do paper
2-SCN-NQ. Cada replica reusa o MESMO pipeline validado (md_simulation.py
--fresh --mode subversive), mas com TAG NOVO -> as velocidades sao reatribuidas
aleatoriamente no warm-up (setVelocitiesToTemperature sem seed) e o integrador
Langevin Middle e estocastico => replicas independentes da MESMA pose inicial.

Ordem (round-robin, menor->maior dentro de cada rodada) para travar o resultado
essencial (n=2 nos 3 sistemas) o quanto antes (~36 h), deixando a 2a rodada
(n=3) como bonus:
    Rodada 1:  hGR_r2, LiTR_r2, TcTR_r2        (~36 h => n=2 completo)
    Rodada 2:  hGR_r3, LiTR_r3, TcTR_r3        (~72 h => n=3 completo)

Robustez
--------
- Subprocesso por corrida: um segfault do OpenMM mata so o filho; o driver segue.
- RESUME automatico: se uma corrida cai com .chk + minimized_<tag>.pdb gravados,
  a retentativa entra em modo --fast --complex minimized_<tag>.pdb, e md_simulation
  retoma a producao do checkpoint (mesmo nº de atomos -> casa). Se nao ha checkpoint
  ainda, retenta --fresh do zero.
- Estado idempotente em replica_state.json: re-lancar o driver pula o que ja terminou.
- Log mestre (replica_driver.log) + log por replica (md_<tag>_run.log).
- Apos cada MD: analyze_binding_mode.py para extrair as metricas da pose
  (distancia ligante-FAD, persistencia) que o paper reporta. Nao-fatal.

NAO chamar com conda activate necessario: usa o python.exe do env datascience
direto e injeta CONDA_PREFIX/PATH no ambiente do subprocesso.
"""
import json
import os
import subprocess
import sys
import time
import datetime
from pathlib import Path

# ----------------------------------------------------------------------------
PY    = r"C:\Users\luizd\miniforge3\envs\datascience\python.exe"
ENVDIR = r"C:\Users\luizd\miniforge3\envs\datascience"
DOCK  = Path(r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking")
MDSIM = DOCK / "md_simulation.py"
KIT   = Path(r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-SNC-NQ"
             r"\montagem_resultados\fig13_3d_kit")
ABM   = KIT / "analyze_binding_mode.py"
WORK  = DOCK / "docking_project" / "figures" / "md_clean"
STATE = DOCK / "replica_state.json"
STATUS = DOCK / "replica_status.json"
MASTER = DOCK / "replica_driver.log"
AGG   = DOCK / "replica_summary.json"

NS = 50.0
PLATFORM = "OpenCL"
MAX_RETRIES = 2          # tentativas extras por corrida (resume)
ABM_STRIDE = "5"

JOBS = [
    {"target": "hGR",  "tag": "hGR_subversive_r2"},
    {"target": "LiTR", "tag": "LiTR_subversive_r2"},
    {"target": "TcTR", "tag": "TcTR_subversive_r2"},
    {"target": "hGR",  "tag": "hGR_subversive_r3"},
    {"target": "LiTR", "tag": "LiTR_subversive_r3"},
    {"target": "TcTR", "tag": "TcTR_subversive_r3"},
]

# Ambiente para os subprocessos (OpenMM acha plugins/DLLs via env do datascience)
ENV = os.environ.copy()
ENV["CONDA_PREFIX"] = ENVDIR
ENV["PATH"] = os.pathsep.join([
    ENVDIR,
    str(Path(ENVDIR) / "Library" / "bin"),
    str(Path(ENVDIR) / "Library" / "mingw-w64" / "bin"),
    str(Path(ENVDIR) / "Scripts"),
    ENV.get("PATH", ""),
])
# Garante que OpenMM use a GPU (OpenCL) e nao limite threads de CPU sem querer
ENV.setdefault("OPENMM_DEFAULT_PLATFORM", "OpenCL")


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    line = "[%s] %s" % (now(), msg)
    print(line, flush=True)
    try:
        with open(MASTER, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:
        pass


def load_json(p, default):
    try:
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        pass
    return default


def save_json(p, obj):
    try:
        p.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    except Exception as e:
        log("  [WARN] nao salvou %s: %s" % (p.name, e))


def md_done(tag):
    """True se a producao terminou: resumo RMSD existe com t_total ~50 ns."""
    j = WORK / ("md_%s_rmsd_summary.json" % tag)
    if not j.exists():
        return False
    try:
        d = json.loads(j.read_text(encoding="utf-8"))
        return float(d.get("t_total_ns", 0)) >= 49.0
    except Exception:
        return False


def build_cmd(tag, target):
    """Escolhe --fresh (1a vez) ou --fast+resume (se ha checkpoint do mesmo sistema)."""
    chk = WORK / ("md_%s.chk" % tag)
    minpdb = WORK / ("minimized_%s.pdb" % tag)
    base = [PY, str(MDSIM),
            "--target", target, "--mode", "subversive",
            "--ns", str(NS), "--tag", tag,
            "--platform", PLATFORM, "--allow-fallback"]
    if chk.exists() and minpdb.exists():
        # retoma do checkpoint usando exatamente o sistema solvatado minimizado
        # (mesmo nº de atomos -> loadCheckpoint casa). md_simulation faz resume.
        return base + ["--fast", "--complex", str(minpdb)], "resume(--fast)"
    return base + ["--fresh"], "fresh"


def run_md_job(tag, target):
    cmd, mode = build_cmd(tag, target)
    runlog = WORK / ("md_%s_run.log" % tag)
    log("  [%s] modo=%s :: %s" % (tag, mode, " ".join(cmd)))
    t0 = time.time()
    try:
        with open(runlog, "a", encoding="utf-8") as fh:
            fh.write("\n\n==== %s :: %s ====\n" % (now(), mode))
            fh.flush()
            rc = subprocess.run(cmd, cwd=str(DOCK), env=ENV,
                                stdout=fh, stderr=subprocess.STDOUT).returncode
    except Exception as e:
        log("  [%s] EXCECAO no subprocess: %s" % (tag, e))
        rc = -99
    dt = (time.time() - t0) / 3600.0
    return rc, dt


def run_binding_mode(tag):
    pdb = WORK / ("minimized_%s.pdb" % tag)
    dcd = WORK / ("md_%s.dcd" % tag)
    if not (pdb.exists() and dcd.exists() and ABM.exists()):
        log("  [%s] binding_mode pulado (inputs/script ausentes)." % tag)
        return
    runlog = WORK / ("binding_mode_%s_run.log" % tag)
    cmd = [PY, str(ABM), str(pdb), str(dcd), tag, ABM_STRIDE]
    log("  [%s] binding_mode :: %s" % (tag, " ".join(cmd)))
    try:
        with open(runlog, "w", encoding="utf-8") as fh:
            subprocess.run(cmd, cwd=str(KIT), env=ENV, timeout=7200,
                           stdout=fh, stderr=subprocess.STDOUT)
    except Exception as e:
        log("  [%s] binding_mode falhou (%s); seguindo." % (tag, e))


def write_status(state, current, t_start):
    elapsed = (time.time() - t_start) / 3600.0
    done = [t for t, s in state.items() if s.get("status") == "done"]
    status = {
        "updated": now(),
        "current_job": current,
        "elapsed_h": round(elapsed, 2),
        "done": done,
        "n_done": len(done),
        "n_total": len(JOBS),
        "jobs": state,
    }
    save_json(STATUS, status)


def aggregate():
    """Junta rmsd_summary + binding_mode de cada replica num arquivo so."""
    agg = {"generated": now(), "replicas": {}}
    # inclui tambem as 3 corridas ORIGINAIS p/ comparacao (n=1 -> n=3)
    originals = ["TcTR_subversive_50ns", "LiTR_subversive", "hGR_subversive"]
    tags = originals + [j["tag"] for j in JOBS]
    for tag in tags:
        entry = {}
        rs = WORK / ("md_%s_rmsd_summary.json" % tag)
        if rs.exists():
            try:
                entry["rmsd_summary"] = json.loads(rs.read_text(encoding="utf-8"))
            except Exception:
                pass
        bm = KIT / "redock" / ("binding_mode_%s.json" % tag)
        if not bm.exists():
            bm = WORK / ("binding_mode_%s.json" % tag)
        if bm.exists():
            try:
                entry["binding_mode"] = json.loads(bm.read_text(encoding="utf-8"))
            except Exception:
                pass
        if entry:
            agg["replicas"][tag] = entry
    save_json(AGG, agg)
    log("  agregado -> %s (%d replicas com dados)" % (AGG.name, len(agg["replicas"])))


def main():
    log("=" * 72)
    log("replica_driver INICIADO | %d corridas | plataforma=%s | %g ns cada"
        % (len(JOBS), PLATFORM, NS))
    log("python=%s" % PY)
    log("workdir=%s" % WORK)
    log("ordem: " + " -> ".join(j["tag"] for j in JOBS))
    log("=" * 72)

    if not MDSIM.exists():
        log("FATAL: md_simulation.py nao encontrado em %s" % MDSIM)
        sys.exit(2)

    state = load_json(STATE, {})
    t_start = time.time()

    for i, job in enumerate(JOBS, 1):
        tag, target = job["tag"], job["target"]
        write_status(state, tag, t_start)

        if state.get(tag, {}).get("status") == "done" or md_done(tag):
            log("[%d/%d] %s: ja concluido; pulando." % (i, len(JOBS), tag))
            state[tag] = {"status": "done"}
            save_json(STATE, state)
            continue

        attempts = int(state.get(tag, {}).get("attempts", 0))
        ok = False
        while attempts <= MAX_RETRIES:
            attempts += 1
            log("[%d/%d] %s: tentativa %d/%d" % (i, len(JOBS), tag, attempts, MAX_RETRIES + 1))
            state[tag] = {"status": "running", "attempts": attempts, "started": now()}
            save_json(STATE, state)
            write_status(state, tag, t_start)

            rc, dt = run_md_job(tag, target)
            done = md_done(tag)
            log("[%d/%d] %s: rc=%s  md_done=%s  (%.1f h)" % (i, len(JOBS), tag, rc, done, dt))

            if rc == 0 and done:
                log("[%d/%d] %s: MD OK. Rodando binding_mode..." % (i, len(JOBS), tag))
                run_binding_mode(tag)
                state[tag] = {"status": "done", "hours": round(dt, 2),
                              "attempts": attempts, "finished": now()}
                save_json(STATE, state)
                ok = True
                break
            else:
                last = attempts > MAX_RETRIES
                state[tag] = {"status": ("failed" if last else "retry"),
                              "attempts": attempts, "rc": rc, "last": now()}
                save_json(STATE, state)
                log("[%d/%d] %s: %s" % (i, len(JOBS), tag,
                    "DESISTINDO deste tag (segue p/ o proximo)" if last
                    else "vai retentar (resume se houver checkpoint)"))

        write_status(state, tag, t_start)
        elapsed = (time.time() - t_start) / 3600.0
        log("  >>> %s: %s | decorrido total: %.1f h"
            % (tag, "OK" if ok else "FALHOU", elapsed))
        # agrega parcialmente a cada corrida (deixa um overview sempre fresco)
        try:
            aggregate()
        except Exception as e:
            log("  [WARN] aggregate parcial falhou: %s" % e)

    log("=" * 72)
    log("TODAS AS CORRIDAS PROCESSADAS. Resumo final:")
    for j in JOBS:
        t = j["tag"]
        log("   %-22s %s" % (t, state.get(t, {}).get("status", "?")))
    aggregate()
    write_status(state, None, t_start)
    log("FIM replica_driver. Total: %.1f h" % ((time.time() - t_start) / 3600.0))


if __name__ == "__main__":
    main()
