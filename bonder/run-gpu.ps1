# =============================================================================
# bonder/run-gpu.ps1 - sobe a plataforma com a GPU (RTX 3060).
# Requer NVIDIA Container Toolkit no WSL2. Uso:
#   powershell -ExecutionPolicy Bypass -File .\run-gpu.ps1
# Depois abra http://localhost:8080
# =============================================================================
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
New-Item -ItemType Directory -Force -Path inputs | Out-Null

Write-Host "== bonder (GPU) ==" -ForegroundColor Cyan
Write-Host "Coloque suas saidas de busca em: $PSScriptRoot\inputs" -ForegroundColor Yellow
Write-Host "Dashboard: http://localhost:8080  (Ctrl+C para parar)" -ForegroundColor Green

docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build
