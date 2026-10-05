# -*- coding: utf-8 -*-
"""
run_md_chain.py -- Orquestrador SERIAL das DMs restantes (1 GPU). Espera a
hGR-subversiva (em curso) terminar e entao roda, em sequencia, cada DM de 50 ns
seguida da analise de modo de ligacao:

  0. (espera) hGR_subversive  -> analyze (ditiol 58,63)
  1. hGR_inhibitor   (mode inhibitor)   -> analyze (58,63)
  2. LiTR_subversive (mode subversive)  -> analyze (52,57)
  3. LiTR_inhibitor  (mode inhibitor)   -> analyze (52,57)

Completa a matriz 2x2 (hGR) + 2 (LiTR), somando aos 2 ja prontos da TcTR.
Sementes/poses devem existir (setup_*_md.py / build_*_inhibitor_pose.py).

Uso (env datascience, background):  python run_md_chain.py
"""
import os, sys, time, subprocess, datetime, json

PY = r"C:\Users\luizd\miniforge3\envs\datascience\python.exe"
DOCK = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
MDSIM = os.path.join(DOCK, "md_simulation.py")
MD = os.path.join(DOCK, "docking_project", "figures", "md_clean")
KIT = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-snc-nq2\montagem_resultados\fig13_3d_kit"
ANALYZE = os.path.join(KIT, "analyze_binding_mode.py")
LOGF = os.path.join(DOCK, "run_md_chain.log")
DONE = os.path.join(DOCK, "run_md_chain.DONE")

# (target, mode, tag, dithiol)
JOBS = [
    ("hGR",  "inhibitor",  "hGR_inhibitor",   "58,63"),
    ("LiTR", "subversive", "LiTR_subversive", "52,57"),
    ("LiTR", "inhibitor",  "LiTR_inhibitor",  "52,57"),
]


def log(m):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), m)
    print(line, flush=True)
    open(LOGF, "a", encoding="utf-8").write(line + "\n")


def summary(tag):
    return os.path.join(MD, "md_%s_rmsd_summary.json" % tag)


def wait_done(tag, poll=120, maxh=20):
    log("aguardando termino de %s ..." % tag)
    t0 = time.time()
    while not os.path.exists(summary(tag)):
        if time.time() - t0 > maxh * 3600:
            log("TIMEOUT esperando %s" % tag); return False
        time.sleep(poll)
    time.sleep(20)
    log("%s concluido." % tag); return True


def analyze(tag, dithiol):
    topo = os.path.join(MD, "minimized_%s.pdb" % tag)
    dcd = os.path.join(MD, "md_%s.dcd" % tag)
    if not (os.path.exists(topo) and os.path.exists(dcd)):
        log("  [analyze %s] topo/dcd ausentes; pulando" % tag); return
    r = subprocess.run([PY, ANALYZE, topo, dcd, tag, "5", dithiol],
                       cwd=KIT, capture_output=True, text=True, encoding="utf-8", errors="replace")
    tail = ((r.stdout or "")[-300:]).replace("\n", " | ")
    log("  [analyze %s] rc=%d  %s" % (tag, r.returncode, tail))


def run_md(target, mode, tag):
    log("=== iniciando DM %s (%s, 50 ns) ===" % (tag, mode))
    cmd = [PY, MDSIM, "--target", target, "--mode", mode, "--fresh",
           "--ns", "50", "--tag", tag, "--platform", "CUDA", "--no-resume"]
    r = subprocess.run(cmd, cwd=DOCK, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    open(os.path.join(DOCK, "md_%s.runlog" % tag), "w", encoding="utf-8").write(
        (r.stdout or "") + "\n" + (r.stderr or ""))
    ok = os.path.exists(summary(tag))
    log("=== DM %s terminou rc=%d (rmsd_summary=%s) ===" % (tag, r.returncode, ok))
    return ok


def main():
    log("==== run_md_chain iniciado ====")
    # 0) hGR_subversive (em curso) -> analise
    if wait_done("hGR_subversive"):
        analyze("hGR_subversive", "58,63")
    status = {"hGR_subversive": os.path.exists(summary("hGR_subversive"))}
    # 1..3) jobs em serie
    for target, mode, tag, dith in JOBS:
        ok = run_md(target, mode, tag)
        status[tag] = ok
        if ok:
            analyze(tag, dith)
        else:
            log("  AVISO: %s nao gerou rmsd_summary; seguindo para o proximo" % tag)
    json.dump({"finished_at": datetime.datetime.now().isoformat(timespec="seconds"),
               "status": status}, open(DONE, "w"), indent=2)
    log("==== CHAIN CONCLUIDA. status=%s ====" % status)
    log("Proximo (Claude): comparar os 6 modos (TcTR/hGR/LiTR x subversivo/inibidor) "
        "-> veredito de seletividade dinamica + figura/tabela do paper.")


if __name__ == "__main__":
    main()
