@echo off
REM ============================================================
REM  md_simulation.bat — Lancador conveniente do md_simulation.py
REM  (antigo md_clean.bat / md_clean.py — renomeado em 2026-05-25)
REM  Usa env conda 'datascience' (miniforge3)
REM  Uso: md_simulation.bat [ns]   (default 5 ns)
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

echo === md_simulation.py --fast --ns %NS% ===
python md_simulation.py --fast --ns %NS%
endlocal
