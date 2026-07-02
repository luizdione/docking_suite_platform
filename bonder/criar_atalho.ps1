# =============================================================================
# criar_atalho.ps1 - cria um atalho "bonder" na Area de Trabalho que inicia o
# app e abre http://localhost:8080. Rode uma vez:
#   powershell -ExecutionPolicy Bypass -File .\criar_atalho.ps1
# =============================================================================
$ErrorActionPreference = "Stop"
$launcher = Join-Path $PSScriptRoot "iniciar_bonder.ps1"
$desktop = [Environment]::GetFolderPath("Desktop")   # respeita a redirecao do OneDrive
$lnk = Join-Path $desktop "bonder.lnk"

$ws = New-Object -ComObject WScript.Shell
$s = $ws.CreateShortcut($lnk)
$s.TargetPath = "powershell.exe"
$s.Arguments = "-ExecutionPolicy Bypass -File `"$launcher`""
$s.WorkingDirectory = $PSScriptRoot
$s.IconLocation = "shell32.dll,13"
$s.Description = "bonder - dashboard CoIP-MS (localhost:8080)"
$s.Save()
Write-Host "Atalho 'bonder' criado em: $lnk" -ForegroundColor Green
