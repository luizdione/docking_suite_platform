# =============================================================================
# bonder/stop.ps1 - derruba os conteineres da plataforma.
#   powershell -ExecutionPolicy Bypass -File .\stop.ps1
# =============================================================================
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
docker compose -f docker/docker-compose.yml down
Write-Host "bonder parado." -ForegroundColor Green
