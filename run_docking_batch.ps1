# run_docking_batch.ps1
# Versao PowerShell do batch de docking — usa env 'datascience' explicitamente
# Uso:
#   .\run_docking_batch.ps1                # docking dos 9 alvos
#   .\run_docking_batch.ps1 -Targets TcTR  # apenas TcTR

param(
    [string]$Targets = "TcTR,LiTR,TcCruzain,TcFeSOD,LmDHODH,hGR,hCatL,hDHODH,hMnSOD",
    [int]$Exhaustiveness = 32
)

$ErrorActionPreference = "Stop"

# Caminho do Python no env 'datascience' (que tem openbabel, rdkit, meeko, vina)
$PythonExe = "C:\Users\luizd\miniforge3\envs\datascience\python.exe"

if (-not (Test-Path $PythonExe)) {
    Write-Host "[ERRO] Env 'datascience' nao encontrado em $PythonExe" -ForegroundColor Red
    Write-Host "Crie o env primeiro com: setup_bioinfo_env.bat"
    exit 1
}

Set-Location $PSScriptRoot

Write-Host ""
Write-Host "========================================================================" -ForegroundColor Green
Write-Host "  DOCKING BATCH - 2-SCN-NQ contra alvos: $Targets" -ForegroundColor Green
Write-Host "========================================================================" -ForegroundColor Green
Write-Host ""
Write-Host "Python EXE: $PythonExe"
Write-Host "Exhaustiveness: $Exhaustiveness"
Write-Host ""

# Verificar dependencias
Write-Host "[1/3] Verificando dependencias..." -ForegroundColor Cyan
& $PythonExe -c "import rdkit; import openbabel; import Bio; print('  [OK] rdkit, openbabel, biopython presentes')" 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "[AVISO] Algumas dependencias podem estar faltando" -ForegroundColor Yellow
    Read-Host "Pressione ENTER para continuar mesmo assim, ou Ctrl+C para abortar"
}

# Rodar docking
Write-Host ""
Write-Host "[2/3] Executando docking (exhaustiveness=$Exhaustiveness)..." -ForegroundColor Cyan
Write-Host "Tempo estimado: ~5-10 min por alvo"
Write-Host ""

& $PythonExe docking_v2.py --targets $Targets --exhaustiveness $Exhaustiveness

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "[ERRO] Docking falhou - veja log em docking_project\logs\" -ForegroundColor Red
    exit 1
}

# Verificar resultados
Write-Host ""
Write-Host "[3/3] Verificando outputs do docking..." -ForegroundColor Cyan
$TargetArray = $Targets -split ","
& $PythonExe check_docking_status.py --targets @TargetArray

Write-Host ""
Write-Host "========================================================================" -ForegroundColor Green
Write-Host "  Docking concluido. Proximo passo: rodar MD com poses docked." -ForegroundColor Green
Write-Host "  Comando sugerido (teste curto 10 ns):" -ForegroundColor Green
Write-Host "    & '$PythonExe' md_module.py --targets TcTR --ns 10 --mode real"
Write-Host "========================================================================" -ForegroundColor Green
