@echo off
REM ============================================================
REM  md_clean.bat — Lancador conveniente do md_clean.py
REM  Usa env conda 'datascience' (miniforge3)
REM  Uso: md_clean.bat [ns]   (default 5 ns)
REM ============================================================
setlocal
cd /d "%~dp0"

set NS=%1
if "%NS%"=="" set NS=5

call "%USERPROFILE%\miniforge3\Scripts\activate.bat" datascience
if errorlevel 1 (
    echo [ERRO] Nao consegui ativar env 'datascience'.
    exit /b 1
)

echo === md_simulation.py --fast --ns %NS% (md_clean.bat e um atalho compat) ===
python md_simulation.py --fast --ns %NS%
endlocal
