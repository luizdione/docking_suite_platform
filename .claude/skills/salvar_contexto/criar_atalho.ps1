<#
  criar_atalho.ps1  -  Cria um atalho (.lnk) do Windows para a pasta desta skill
  (salvar_contexto), para acesso rapido a partir da sua pasta de Projetos.

  O ALVO do atalho e a propria pasta onde este script esta (a pasta da skill),
  resolvida por $PSScriptRoot. Assim funciona independentemente do nome da pasta
  local do repo (docking / docking_suite).

  Uso (no Windows, PowerShell):
    .\criar_atalho.ps1                          # cria em C:\Users\luizd\OneDrive\Data_Science\Projetos
    .\criar_atalho.ps1 -DestDir "D:\OutroLugar" # muda a pasta de destino
    .\criar_atalho.ps1 -Nome "skill_memoria"    # muda o nome do atalho

  NAO requer Administrador.
#>
[CmdletBinding()]
param(
  [string]$DestDir = "C:\Users\luizd\OneDrive\Data_Science\Projetos",
  [string]$Nome    = "salvar_contexto (skill)"
)

$ErrorActionPreference = "Stop"

# Alvo do atalho = a pasta desta skill (onde este script esta).
$Alvo = $PSScriptRoot

if (-not (Test-Path -LiteralPath $DestDir)) {
  Write-Host " Pasta de destino nao existe:" -ForegroundColor Red
  Write-Host "   $DestDir" -ForegroundColor Red
  Write-Host " Ajuste o parametro -DestDir e rode de novo." -ForegroundColor Yellow
  exit 1
}

$Lnk = Join-Path $DestDir ($Nome + ".lnk")

$ws = New-Object -ComObject WScript.Shell
$sc = $ws.CreateShortcut($Lnk)
$sc.TargetPath       = $Alvo
$sc.WorkingDirectory = $Alvo
$sc.Description       = "Skill /salvar_contexto - memoria de sessao e handoff"
$sc.Save()

Write-Host "===================================================================" -ForegroundColor Cyan
Write-Host " Atalho da skill /salvar_contexto criado" -ForegroundColor Green
Write-Host "===================================================================" -ForegroundColor Cyan
Write-Host "   Atalho : $Lnk" -ForegroundColor Green
Write-Host "   Aponta : $Alvo" -ForegroundColor DarkGray
Write-Host ""
Write-Host " Abra a pasta Projetos e clique no atalho para ir direto a skill." -ForegroundColor Yellow
