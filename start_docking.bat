@echo off
title NQ-Docking Suite Dashboard
REM =============================================================================
REM Inicia o Dashboard usando o ambiente conda 'datascience' do Miniforge3
REM Conda: C:\Users\luizd\miniforge3\envs\datascience
REM =============================================================================

cd /d C:\Users\luizd\OneDrive\Data_Science\Projetos\docking

REM Ativa o conda base e depois o ambiente datascience
call C:\Users\luizd\miniforge3\Scripts\activate.bat C:\Users\luizd\miniforge3
call conda activate datascience

if %ERRORLEVEL% neq 0 (
    echo [ERRO] Ambiente 'datascience' nao encontrado!
    echo Execute setup_bioinfo_env.bat primeiro.
    pause
    exit /b 1
)

echo.
echo Ambiente: datascience (Miniforge3)
python --version
echo.

python dashboard_v2.py
