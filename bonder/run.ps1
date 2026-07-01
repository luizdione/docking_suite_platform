# =============================================================================
# bonder/run.ps1 - sobe a plataforma em CPU (funciona em qualquer maquina).
# Uso (PowerShell, dentro da pasta bonder):
#   powershell -ExecutionPolicy Bypass -File .\run.ps1
# Depois abra http://localhost:8080
# =============================================================================
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
New-Item -ItemType Directory -Force -Path inputs | Out-Null

Write-Host "== bonder ==" -ForegroundColor Cyan
Write-Host "Coloque suas saidas de busca em: $PSScriptRoot\inputs" -ForegroundColor Yellow
Write-Host "Subindo (CPU). A primeira vez baixa/compila as imagens (alguns minutos)..." -ForegroundColor Cyan
Write-Host "Dashboard: http://localhost:8080  (Ctrl+C para parar)" -ForegroundColor Green

docker compose -f docker/docker-compose.yml up --build
