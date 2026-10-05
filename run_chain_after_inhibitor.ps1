<#
run_chain_after_inhibitor.ps1
Encadeia automaticamente, SEM supervisao:
  1) Aguarda o MD inibidor 50ns terminar
  2) Roda o post_50ns_pipeline.py (analises + figuras + build paper)
  3) Inicia o MD subversivo 50ns fresh (face-si do FAD, tag TcTR_subversive_50ns)

Rodar como ADMIN (para o power-cap):
  Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
  cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
  .\run_chain_after_inhibitor.ps1
#>

$DOCK   = "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
$Py     = "$env:USERPROFILE\miniforge3\envs\datascience\python.exe"
$SUMMARY= "$DOCK\docking_project\figures\md_clean\md_TcTR_inhibitor_rmsd_summary.json"
$PIPE   = "C:\Users\luizd\OneDrive\Data_Science\Projetos\2-snc-nq2\montagem_resultados\fig13_3d_kit\post_50ns_pipeline.py"
$LIG_SUB= "$DOCK\docking_project\figures\smoke_test_fad\lig_TcTR_fadsite.pdb"
$LOG    = "$DOCK\chain_log.txt"

function Log($msg) {
    $line = "[$((Get-Date).ToString('HH:mm:ss'))] $msg"
    Write-Host $line
    Add-Content $LOG $line
}

Log "=== CHAIN INICIADO ==="
Log "Aguardando fim do inibidor 50ns ($SUMMARY) ..."

# 1) Aguarda o inibidor terminar (poll a cada 3 min, timeout 12h)
$waited = 0
while (-not (Test-Path $SUMMARY)) {
    Start-Sleep 180
    $waited += 3
    if ($waited % 30 -eq 0) {
        $last = (Get-Content "$DOCK\docking_project\figures\md_clean\md_TcTR_inhibitor.csv" -Tail 1) -split ","
        Log "  inibidor ainda rodando: progresso=$($last[0]) speed=$($last[9]) ns/dia"
    }
    if ($waited -ge 720) { Log "TIMEOUT 12h sem conclusao"; exit 1 }
}
Log "Inibidor 50ns CONCLUIDO."

# 2) Analises + paper build
Log "Rodando post_50ns_pipeline.py ..."
& $Py $PIPE
Log "Pipeline concluido (rc=$LASTEXITCODE)."

# 3) MD subversivo 50ns fresh (tag novo para nao sobrescrever o 20ns)
Log "Iniciando MD subversivo 50ns (TcTR_subversive_50ns) ..."
nvidia-smi -pl 130 | Out-Null
Set-Location $DOCK
& $Py md_simulation.py `
    --target TcTR --ns 50 --tag TcTR_subversive_50ns --mode subversive `
    --platform CUDA --allow-fallback --no-resume --fresh `
    --lig-pdb $LIG_SUB

Log "Subversivo 50ns CONCLUIDO (rc=$LASTEXITCODE)."
Log "=== CHAIN ENCERRADO ==="
