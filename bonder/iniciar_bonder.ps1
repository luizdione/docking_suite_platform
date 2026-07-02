# =============================================================================
# iniciar_bonder.ps1 - inicia o dashboard bonder e abre o navegador.
# Usa o Python do conda 'datascience'. Uma janela do Python (servidor) fica
# aberta; feche-a para encerrar o app. Alvo do atalho da Area de Trabalho.
# =============================================================================
$ErrorActionPreference = "SilentlyContinue"
Set-Location -Path $PSScriptRoot

$py = "C:\Users\luizd\miniforge3\envs\datascience\python.exe"
if (-not (Test-Path $py)) { $py = "python" }

# sobe o servidor em uma janela propria (nao bloqueia este script)
Start-Process -FilePath $py -ArgumentList "-m", "bonder.web.app" -WorkingDirectory $PSScriptRoot
# da tempo do servidor subir e abre o navegador
Start-Sleep -Seconds 4
Start-Process "http://localhost:8080"
