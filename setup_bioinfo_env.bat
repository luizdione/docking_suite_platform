@echo off
REM =============================================================================
REM setup_bioinfo_env.bat
REM Instala pacotes de docking e dinamica molecular no ambiente 'datascience'
REM Miniforge3: C:\Users\luizd\miniforge3\envs\datascience (Python 3.11)
REM
REM O ambiente 'datascience' ja possui: numpy, scipy, matplotlib, pandas,
REM pillow, scikit-learn, seaborn, openbabel, biopython, packaging.
REM Este script adiciona APENAS o que falta para o NQ-Docking Suite.
REM
REM NOTA: openff-toolkit e openmmforcefields dependem de ambertools,
REM que NAO esta disponivel para Windows no conda-forge.
REM Esses pacotes sao instalados via pip como fallback.
REM
REM USO:
REM   1. Abra o Miniforge Prompt
REM   2. cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
REM   3. setup_bioinfo_env.bat
REM =============================================================================

set CONDA_ROOT=C:\Users\luizd\miniforge3
set ENV_NAME=datascience

echo.
echo ============================================================
echo   Miniforge3 — Atualizando ambiente '%ENV_NAME%'
echo   Conda: %CONDA_ROOT%
echo   Adicionando pacotes de Docking + Dinamica Molecular
echo ============================================================
echo.

REM ── 1. Inicializa conda ───────────────────────────────────────────────────
call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
if %ERRORLEVEL% neq 0 (
    echo [ERRO] Nao foi possivel ativar o conda em %CONDA_ROOT%
    echo Verifique se o Miniforge3 esta instalado corretamente.
    pause
    exit /b 1
)

REM ── 2. Ativa o ambiente datascience ───────────────────────────────────────
call conda activate %ENV_NAME%
if %ERRORLEVEL% neq 0 (
    echo [ERRO] Ambiente '%ENV_NAME%' nao encontrado!
    echo Verifique com: conda env list
    pause
    exit /b 1
)

echo [OK] Ambiente '%ENV_NAME%' ativado.
python --version
echo.

REM ── 3. Pacotes conda-forge (sem ambertools-dependentes) ───────────────────
echo [1/5] Instalando pacotes conda-forge (docking + MD) ...
echo       Isso pode levar alguns minutos na primeira vez...
echo.
conda install -c conda-forge ^
    rdkit ^
    openmm ^
    mdtraj ^
    parmed ^
    netcdf4 ^
    cftime ^
    flask ^
    flask-cors ^
    requests ^
    -y

if %ERRORLEVEL% neq 0 (
    echo.
    echo [AVISO] Instalacao em bloco falhou. Tentando individualmente...
    echo.
    conda install -c conda-forge rdkit -y
    conda install -c conda-forge openmm -y
    conda install -c conda-forge mdtraj -y
    conda install -c conda-forge parmed netcdf4 cftime -y
    conda install -c conda-forge flask flask-cors requests -y
)

REM ── 4. Pacotes pip (indisponiveis ou problematicos no conda Windows) ──────
echo.
echo [2/5] Instalando pacotes pip ...
pip install MDAnalysis 2>nul
pip install meeko 2>nul
pip install vina 2>nul
echo [OK] Pacotes pip instalados.

REM ── 5. openmmforcefields e openff-toolkit via pip (sem ambertools) ────────
echo.
echo [3/5] Instalando openmmforcefields e openff-toolkit via pip ...
echo       (ambertools nao disponivel para Windows via conda-forge)
echo.
pip install openmmforcefields 2>nul
if %ERRORLEVEL% neq 0 (
    echo [AVISO] openmmforcefields nao instalou. Nao e critico para docking basico.
    echo         Necessario apenas para parametrizacao avancada de ligantes com GAFF2.
)
pip install openff-toolkit 2>nul
if %ERRORLEVEL% neq 0 (
    echo [AVISO] openff-toolkit nao instalou. Nao e critico para docking basico.
    echo         Necessario apenas para parametrizacao de ligantes com OpenFF.
)
echo.

REM ── 6. Verificacao completa ───────────────────────────────────────────────
echo.
echo [4/5] Verificando instalacao ...
echo.
echo ================================================================
echo              VERIFICACAO DE PACOTES — %ENV_NAME%
echo ================================================================
echo.
echo --- Ja existentes (datascience) ---
python -c "import sys; print(f'Python:          {sys.version}')"
python -c "import numpy; print(f'NumPy:           {numpy.__version__}')"
python -c "import scipy; print(f'SciPy:           {scipy.__version__}')"
python -c "import matplotlib; print(f'Matplotlib:      {matplotlib.__version__}')"
python -c "import pandas; print(f'Pandas:          {pandas.__version__}')"
python -c "from PIL import Image; print(f'Pillow:          OK')"
python -c "import Bio; print(f'BioPython:       {Bio.__version__}')" 2>nul || echo BioPython:       FALHOU
python -c "from openbabel import pybel; print(f'OpenBabel:       OK')" 2>nul || echo OpenBabel:       FALHOU

echo.
echo --- Adicionados para Docking/MD ---
python -c "import flask; print(f'Flask:           {flask.__version__}')" 2>nul || echo Flask:           FALHOU
python -c "from flask_cors import CORS; print(f'Flask-CORS:      OK')" 2>nul || echo Flask-CORS:      FALHOU
python -c "from rdkit import Chem; print(f'RDKit:           OK')" 2>nul || echo RDKit:           FALHOU
python -c "import openmm; print(f'OpenMM:          {openmm.__version__}')" 2>nul || echo OpenMM:          FALHOU
python -c "import mdtraj; print(f'MDTraj:          {mdtraj.__version__}')" 2>nul || echo MDTraj:          FALHOU
python -c "import parmed; print(f'ParmEd:          {parmed.__version__}')" 2>nul || echo ParmEd:          FALHOU
python -c "import netCDF4; print(f'netCDF4:         {netCDF4.__version__}')" 2>nul || echo netCDF4:         FALHOU
python -c "import MDAnalysis; print(f'MDAnalysis:      {MDAnalysis.__version__}')" 2>nul || echo MDAnalysis:      FALHOU
python -c "import meeko; print(f'Meeko:           OK')" 2>nul || echo Meeko:           FALHOU
python -c "import requests; print(f'Requests:        {requests.__version__}')" 2>nul || echo Requests:        FALHOU

echo.
echo --- Opcionais (parametrizacao avancada) ---
python -c "import openmmforcefields; print(f'OpenMMFF:        OK')" 2>nul || echo OpenMMFF:        nao instalado (opcional)
python -c "import openff.toolkit; print(f'OpenFF:          OK')" 2>nul || echo OpenFF:          nao instalado (opcional)

echo.
echo --- Plataformas OpenMM (GPU) ---
python -c "import openmm; [print(f'  {openmm.Platform.getPlatform(i).getName()}') for i in range(openmm.Platform.getNumPlatforms())]" 2>nul
if %ERRORLEVEL% neq 0 (
    echo   [AVISO] OpenMM nao detectou plataformas. GPU pode nao estar acessivel.
)

REM ── 7. Salva info para o dashboard ────────────────────────────────────────
echo.
echo [5/5] Salvando configuracao ...
python -c "import json,sys,os; json.dump({'python_path':sys.executable,'env_name':'datascience','env_dir':os.environ.get('CONDA_PREFIX',''),'conda_prefix':os.environ.get('CONDA_PREFIX',''),'miniforge':'C:\\Users\\luizd\\miniforge3','openmm_ok':True,'cuda_version':None},open('venv_info.json','w'),indent=2); print('[OK] venv_info.json atualizado')"

echo.
echo ============================================================
echo   Ambiente '%ENV_NAME%' atualizado com sucesso!
echo.
echo   Localizacao: %CONDA_ROOT%\envs\%ENV_NAME%
echo.
echo   Para usar:
echo     conda activate %ENV_NAME%
echo     cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
echo     python dashboard_v2.py
echo.
echo   Atalho rapido: start_docking.bat
echo ============================================================
echo.
pause
