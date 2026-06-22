#!/usr/bin/env python3
"""
=============================================================================
setup_openmm_venv.py
Instala/verifica OpenMM + dependências de MD no ambiente conda 'datascience'.
Miniforge3: C:\\Users\\luizd\\miniforge3
RTX 3060 / CUDA 12.x — AMD Ryzen 5 5500 / 32 GB RAM

Uso:
    # Ative o ambiente bioinfo primeiro:
    conda activate datascience

    # Depois execute:
    python setup_openmm_venv.py
=============================================================================
"""
import subprocess
import sys
import os
import json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

# ── Detecta o ambiente Python atual (conda ou venv) ─────────────────────────
CONDA_PREFIX = os.environ.get("CONDA_PREFIX", "")
MINIFORGE_ROOT = Path(r"C:\Users\luizd\miniforge3")

if CONDA_PREFIX:
    # Rodando dentro de um ambiente conda
    PYTHON = sys.executable
    PIP = [sys.executable, "-m", "pip"]
    ENV_NAME = os.environ.get("CONDA_DEFAULT_ENV", "unknown")
    ENV_DIR = Path(CONDA_PREFIX)
    print(f"Ambiente conda: {ENV_NAME}")
    print(f"Prefix:         {CONDA_PREFIX}")
else:
    # Fallback: tenta localizar o ambiente bioinfo
    bioinfo_dir = MINIFORGE_ROOT / "envs" / "datascience"
    if sys.platform == "win32":
        bioinfo_py = bioinfo_dir / "python.exe"
    else:
        bioinfo_py = bioinfo_dir / "bin" / "python"

    if bioinfo_py.exists():
        PYTHON = str(bioinfo_py)
        PIP = [str(bioinfo_py), "-m", "pip"]
        ENV_NAME = "datascience"
        ENV_DIR = bioinfo_dir
        print(f"Usando datascience encontrado em: {bioinfo_dir}")
    else:
        # Ultimo fallback: Python atual
        PYTHON = sys.executable
        PIP = [sys.executable, "-m", "pip"]
        ENV_NAME = "sistema"
        ENV_DIR = Path(sys.prefix)
        print(f"[AVISO] Ambiente 'datascience' nao encontrado.")
        print(f"        Usando Python do sistema: {PYTHON}")
        print(f"        Recomendado: execute setup_bioinfo_env.bat primeiro.")

print(f"Python: {PYTHON}")
print()


def run(cmd: list, check: bool = True) -> int:
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    if check and result.returncode != 0:
        print(f"  [ERRO] Código {result.returncode}")
    return result.returncode


def detect_cuda() -> str | None:
    """Detecta versão do CUDA disponível."""
    try:
        r = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
        if r.returncode == 0:
            import re
            m = re.search(r"CUDA Version:\s*([\d.]+)", r.stdout)
            if m:
                return m.group(1)
        return None
    except FileNotFoundError:
        return None


def check_openmm() -> bool:
    """Verifica se OpenMM já está instalado."""
    try:
        r = subprocess.run([PYTHON, "-c", "import openmm; print(openmm.__version__)"],
                           capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout.strip()
        return None
    except Exception:
        return None


# ═════════════════════════════════════════════════════════════════════════════
print("=" * 62)
print("  SETUP OpenMM + MD Stack")
print(f"  Ambiente: {ENV_NAME} (Miniforge3)")
print("  RTX 3060 12 GB | CUDA | AMD Ryzen 5 5500 | 32 GB RAM")
print("=" * 62)
print()

# 1. Verificar CUDA
cuda_ver = detect_cuda()
if cuda_ver:
    print(f"  [OK] CUDA detectado: {cuda_ver}")
    print(f"       RTX 3060 suportada — plataforma CUDA disponivel")
else:
    print("  [AVISO] nvidia-smi nao encontrado. Verifique os drivers CUDA.")
    print("         OpenMM ainda pode instalar mas rodara em CPU.")
print()

# 2. Verificar se já instalado
existing = check_openmm()
if existing:
    print(f"  [OK] OpenMM {existing} ja instalado!")
    print("       Para reinstalar, delete o pacote e re-execute.")
    print()
    REINSTALL = input("  Reinstalar mesmo assim? [s/N]: ").strip().lower()
    if REINSTALL != "s":
        print("  Abortando. OpenMM ja esta pronto.")
        sys.exit(0)
print()

# 3. Detectar se estamos em conda — se sim, usar conda install
is_conda = bool(CONDA_PREFIX)

if is_conda:
    # ── Via conda (metodo preferido para GPU) ──
    print("─── Passo 1: Instalando OpenMM via conda-forge (recomendado) ───")
    run(["conda", "install", "-c", "conda-forge", "openmm", "-y"], check=False)
    print()

    print("─── Passo 2: Instalando MDTraj via conda-forge ───")
    run(["conda", "install", "-c", "conda-forge", "mdtraj", "-y"], check=False)
    print()

    print("─── Passo 3: Instalando OpenFF + OpenMMForceFields via conda-forge ───")
    run(["conda", "install", "-c", "conda-forge",
         "openff-toolkit", "openmmforcefields", "-y"], check=False)
    print()

    print("─── Passo 4: Instalando ParmEd + netCDF4 via conda-forge ───")
    run(["conda", "install", "-c", "conda-forge", "parmed", "netcdf4", "cftime", "-y"],
        check=False)
    print()

    print("─── Passo 5: Verificando RDKit ───")
    rc = run([PYTHON, "-c", "from rdkit import Chem; print('OK')"], check=False)
    if rc != 0:
        run(["conda", "install", "-c", "conda-forge", "rdkit", "-y"], check=False)
    print()

    print("─── Passo 6: Instalando MDAnalysis via pip ───")
    run(PIP + ["install", "MDAnalysis"], check=False)
    print()

else:
    # ── Via pip (fallback se nao esta em conda) ──
    print("─── Passo 1: Atualizando pip + setuptools ───")
    run(PIP + ["install", "--upgrade", "pip", "setuptools", "wheel"])
    print()

    print("─── Passo 2: Instalando OpenMM (via pip) ───")
    rc = run(PIP + ["install", "openmm", "--prefer-binary"], check=False)
    if rc != 0:
        print()
        print("  *** RECOMENDACAO ***")
        print("  Para RTX 3060 + CUDA 12, use o ambiente conda centralizado:")
        print()
        print("    Abra o Miniforge Prompt e execute:")
        print("    cd C:\\Users\\luizd\\OneDrive\\Data_Science\\Projetos\\docking")
        print("    setup_bioinfo_env.bat")
        print()
        sys.exit(1)
    print()

    print("─── Passo 3: Instalando MDTraj ───")
    run(PIP + ["install", "mdtraj", "--prefer-binary"], check=False)
    print()

    print("─── Passo 4: Instalando OpenFF + OpenMMForceFields ───")
    run(PIP + ["install", "openff-toolkit", "openmmforcefields", "--prefer-binary"],
        check=False)
    print()

    print("─── Passo 5: Instalando ParmEd ───")
    run(PIP + ["install", "parmed"], check=False)
    print()

    print("─── Passo 6: Verificando RDKit ───")
    rc = run([PYTHON, "-c", "from rdkit import Chem; print('OK')"], check=False)
    if rc != 0:
        run(PIP + ["install", "rdkit", "--prefer-binary"], check=False)
    print()

# ── Verificação final ──
print("─── Verificação final ───")
tests = [
    ("openmm",     "import openmm; print('OpenMM:', openmm.__version__)"),
    ("cuda",       "import openmm; p=openmm.Platform.getPlatformByName('CUDA'); "
                   "print('CUDA OK:', p.getName())"),
    ("opencl",     "import openmm; p=openmm.Platform.getPlatformByName('OpenCL'); "
                   "print('OpenCL OK:', p.getName())"),
    ("mdtraj",     "import mdtraj; print('MDTraj:', mdtraj.__version__)"),
    ("rdkit",      "from rdkit import Chem; print('RDKit OK')"),
    ("numpy",      "import numpy; print('NumPy:', numpy.__version__)"),
    ("matplotlib", "import matplotlib; print('Matplotlib:', matplotlib.__version__)"),
    ("MDAnalysis", "import MDAnalysis; print('MDAnalysis:', MDAnalysis.__version__)"),
]

all_ok = True
for name, code in tests:
    r = subprocess.run([PYTHON, "-c", code], capture_output=True, text=True)
    if r.returncode == 0:
        print(f"  [OK] {name}: {r.stdout.strip()}")
    else:
        if name in ("cuda", "opencl"):
            print(f"  [--] {name}: nao disponivel (nao critico)")
        else:
            print(f"  [ERRO] {name}: {r.stderr.strip()[:120]}")
            all_ok = False

print()
if all_ok:
    print("=" * 62)
    print("  INSTALACAO COMPLETA!")
    print(f"  Ambiente: {ENV_NAME}")
    print()
    print("  No dashboard, vá em MD Simulation →")
    print("  Modo: 'OpenMM Real' e clique 'Verificar OpenMM'.")
    print("=" * 62)
else:
    print("=" * 62)
    print("  INSTALACAO PARCIAL — veja erros acima.")
    print()
    print("  Se GPU nao disponivel, o dashboard usa CPU como fallback.")
    print("=" * 62)

# Salvar informações para o dashboard
venv_info_path = SCRIPT_DIR / "venv_info.json"
with open(venv_info_path, "w") as f:
    json.dump({
        "python_path":  PYTHON,
        "env_name":     ENV_NAME,
        "env_dir":      str(ENV_DIR),
        "conda_prefix": CONDA_PREFIX,
        "miniforge":    str(MINIFORGE_ROOT),
        "openmm_ok":    all_ok,
        "cuda_version": cuda_ver,
    }, f, indent=2)
print(f"\n  Configuracao salva em: {venv_info_path}")
print("  O dashboard usa esse arquivo para encontrar o Python MD.")
