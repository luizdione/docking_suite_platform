#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gpu_gate.py - PORTAO de GPU. Espera o LILO/Ollama LIBERAR a GPU e SO ENTAO
dispara replica_driver.py (destacado). Evita contender com o LLM e derrubar o
servidor.

Sinal de "GPU livre" (robusto, sem precisar de admin):
  (1) nenhum processo de compute chamado 'llama-server'/'ollama' na GPU  E
  (2) VRAM usada < FREE_MB (garante folga p/ a MD ~2 GB)
sustentado por IDLE_NEEDED checagens consecutivas (evita disparar num intervalo
curto entre inferencias do LILO).

Lanca o replica_driver.py de forma DESTACADA (sobrevive ao fim deste gate) e sai.
Se a GPU nunca liberar dentro de MAX_WAIT_H, NAO dispara (evita OOM) e loga
instrucao para liberar o LILO manualmente (ex.: 'ollama stop <modelo>').
"""
import subprocess
import time
import datetime
import json
import sys
from pathlib import Path

DOCK   = Path(r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking")
PY     = r"C:\Users\luizd\miniforge3\envs\datascience\python.exe"
DRIVER = DOCK / "replica_driver.py"
LOG    = DOCK / "gpu_gate.log"
STATUS = DOCK / "gpu_gate_status.json"
SMI    = r"C:\Windows\System32\nvidia-smi.exe"

FREE_MB     = 7000     # VRAM usada abaixo disso (com LILO fora) = folga p/ MD
POLL_S      = 300      # 5 min entre checagens
IDLE_NEEDED = 3        # nº de checagens livres consecutivas (~15 min)
MAX_WAIT_H  = 70       # teto de espera (dentro da janela de ~3 dias)
LILO_NAMES  = ("llama", "ollama")


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(m):
    line = "[%s] %s" % (now(), m)
    print(line, flush=True)
    try:
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:
        pass


def smi(args):
    exe = SMI if Path(SMI).exists() else "nvidia-smi"
    return subprocess.run([exe] + args, capture_output=True, text=True, timeout=60)


def gpu_state():
    """Retorna (mem_used_mb, util_pct, lilo_present)."""
    try:
        r = smi(["--query-gpu=memory.used,utilization.gpu",
                 "--format=csv,noheader,nounits"])
        mem, util = [int(x.strip()) for x in r.stdout.strip().split(",")[:2]]
    except Exception as e:
        log("  [WARN] query-gpu falhou: %s; assumindo OCUPADA." % e)
        return 999999, 100, True
    try:
        a = smi(["--query-compute-apps=name", "--format=csv,noheader"])
        apps = a.stdout.lower()
        lilo = any(n in apps for n in LILO_NAMES)
    except Exception as e:
        log("  [WARN] query-compute-apps falhou: %s; assumindo LILO presente." % e)
        lilo = True
    return mem, util, lilo


def launch_driver():
    out = open(DOCK / "replica_driver.boot.out.log", "a", encoding="utf-8")
    err = open(DOCK / "replica_driver.boot.err.log", "a", encoding="utf-8")
    DETACHED = 0x00000008
    NEWGROUP = 0x00000200
    p = subprocess.Popen([PY, str(DRIVER)], cwd=str(DOCK),
                         stdout=out, stderr=err,
                         creationflags=DETACHED | NEWGROUP, close_fds=True)
    try:
        (DOCK / "replica_driver.pid").write_text(str(p.pid))
    except Exception:
        pass
    log(">>> GPU LIVRE confirmada. DISPAREI replica_driver.py (PID %d)." % p.pid)
    return p.pid


def write_status(mem, util, lilo, free, idle, el, fired=False, pid=None):
    try:
        STATUS.write_text(json.dumps({
            "updated": now(), "vram_used_mb": mem, "util_pct": util,
            "lilo_present": lilo, "free": free,
            "idle_streak": idle, "need": IDLE_NEEDED,
            "elapsed_h": round(el, 2), "fired": fired, "driver_pid": pid,
        }, indent=2), encoding="utf-8")
    except Exception:
        pass


def main():
    log("=" * 60)
    log("gpu_gate INICIADO. Dispara a MD so quando LILO sair da GPU e "
        "VRAM<%dMiB por %dx%ds." % (FREE_MB, IDLE_NEEDED, POLL_S))
    if not DRIVER.exists():
        log("FATAL: replica_driver.py ausente em %s" % DRIVER)
        sys.exit(2)
    t0 = time.time()
    idle = 0
    while True:
        mem, util, lilo = gpu_state()
        free = (not lilo) and (mem < FREE_MB)
        idle = idle + 1 if free else 0
        el = (time.time() - t0) / 3600.0
        write_status(mem, util, lilo, free, idle, el)
        log("VRAM=%dMiB util=%d%% LILO=%s livre=%s streak=%d/%d t=%.1fh"
            % (mem, util, lilo, free, idle, IDLE_NEEDED, el))
        if idle >= IDLE_NEEDED:
            pid = launch_driver()
            write_status(mem, util, lilo, free, idle, el, fired=True, pid=pid)
            break
        if el > MAX_WAIT_H:
            log("TETO de %dh sem GPU livre. NAO disparei (evita OOM). "
                "Para rodar: libere o LILO (ex.: 'ollama stop') e re-lance "
                "gpu_gate.py, ou rode replica_driver.py manualmente." % MAX_WAIT_H)
            break
        time.sleep(POLL_S)
    log("FIM gpu_gate.")


if __name__ == "__main__":
    main()
