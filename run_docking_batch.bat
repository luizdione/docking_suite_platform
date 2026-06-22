@echo off
REM ============================================================================
REM run_docking_batch.bat - Roda docking molecular para os 9 alvos do estudo
REM ============================================================================
REM Saida: docking_project\docking\result_<TID>.pdbqt para cada alvo
REM Tempo estimado: ~5-10 min por alvo (exhaustiveness=32) = ~1.5 horas total
REM ============================================================================

setlocal

REM Usar SEMPRE o env 'datascience' (tem openbabel, rdkit, meeko, vina, openmm)
REM NAO usar o env 'base' que nao tem essas dependencias.
set PYTHON_EXE=C:\Users\luizd\miniforge3\envs\datascience\python.exe
if not exist "%PYTHON_EXE%" (
    echo [ERRO] Env 'datascience' nao encontrado em C:\Users\luizd\miniforge3\envs\datascience
    echo Crie o env primeiro com: setup_bioinfo_env.bat
    pause
    exit /b 1
)

cd /d "%~dp0"

REM Verificar se as dependencias criticas estao instaladas
echo.
echo Verificando dependencias do env 'datascience'...
"%PYTHON_EXE%" -c "import sys; ok=True; [exec('try:\n    __import__(m)\n    print(f\"  [OK] {m}\")\nexcept ImportError:\n    print(f\"  [FALTA] {m}\"); ok=False') for m in ['rdkit','openbabel','meeko','Bio']]; sys.exit(0 if ok else 1)"
if errorlevel 1 (
    echo.
    echo [AVISO] Algumas dependencias estao faltando no env 'datascience'.
    echo O docking ainda pode funcionar com fallbacks. Continuar mesmo assim?
    pause
)

echo.
echo ========================================================================
echo   DOCKING BATCH - 2-SCN-NQ contra 9 alvos (5 parasita + 4 humanos)
echo ========================================================================
echo.
echo Alvos parasitarios (5):
echo   TcTR       1BZL  Trypanothione Reductase    (T. cruzi)     [FAD]
echo   LiTR       2JK6  Trypanothione Reductase    (L. infantum)  [FAD]
echo   TcCruzain  3KKU  Cruzain                    (T. cruzi)
echo   TcFeSOD    4DVH  Fe-Superoxide Dismutase    (T. cruzi)     [Fe]
echo   LmDHODH    3MHU  Dihydroorotate Dehydrog.   (L. major)     [FMN+FAD]
echo.
echo Controles humanos (4):
echo   hGR        3GRS  Glutathione Reductase      (H. sapiens)   [FAD]
echo   hCatL      2XU3  Cathepsin L                (H. sapiens)
echo   hDHODH     1D3H  Dihydroorotate Dehydrog.   (H. sapiens)   [FMN]
echo   hMnSOD     1LUV  Mn-Superoxide Dismutase    (H. sapiens)   [Mn]
echo.
echo Pre-requisitos:
echo   - vina.exe presente nesta pasta
echo   - PDBs serao baixados automaticamente do RCSB
echo   - Ligante 2-SCN-NQ ja preparado em docking_project\ligand\
echo.
pause

REM Etapa 1: Listar alvos disponiveis (sanity check)
echo.
echo [1/3] Listando alvos suportados...
"%PYTHON_EXE%" docking_v2.py --list
echo.

REM Etapa 2: Rodar docking
echo.
echo [2/3] Executando docking (exhaustiveness=32)...
echo Isso pode demorar ~1.5 horas total.
echo.
"%PYTHON_EXE%" docking_v2.py --targets TcTR,LiTR,TcCruzain,TcFeSOD,LmDHODH,hGR,hCatL,hDHODH,hMnSOD --exhaustiveness 32

if errorlevel 1 (
    echo.
    echo [ERRO] Docking falhou - veja log em docking_project\logs\
    pause
    exit /b 1
)

REM Etapa 3: Verificar resultados
echo.
echo [3/3] Verificando outputs do docking...
"%PYTHON_EXE%" check_docking_status.py --targets TcTR LiTR TcCruzain TcFeSOD LmDHODH hGR hCatL hDHODH hMnSOD

echo.
echo ========================================================================
echo   Docking concluido. Proximo passo: rodar MD com poses docked.
echo   Comando sugerido (teste curto 10 ns):
echo     %PYTHON_EXE% md_module.py --targets TcTR --ns 10 --mode real
echo ========================================================================
pause
endlocal
