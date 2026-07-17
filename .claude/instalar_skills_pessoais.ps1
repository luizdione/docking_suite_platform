<#
  instalar_skills_pessoais.ps1  -  Instala/atualiza as skills pessoais no Windows.

  O que faz:
    * copia as skills salvar_contexto e salvar_contexto_compactar (deste repo)
      para  %USERPROFILE%\.claude\skills\  (escopo pessoal = vale em TODOS os projetos);
    * copia o manual_skills.pdf para a sua pasta Projetos;
    * (opcional -HookGlobal) instala o hook de auto-load do CONTEXTO_SESSAO.md no
      settings.json GLOBAL, para funcionar em qualquer projeto.

  Uso (PowerShell, na pasta .claude do repo):
    .\instalar_skills_pessoais.ps1
    .\instalar_skills_pessoais.ps1 -HookGlobal          # + hook global (todos os projetos)
    .\instalar_skills_pessoais.ps1 -ProjetosDir "D:\..."# muda o destino do PDF

  NAO requer Administrador.
#>
[CmdletBinding()]
param(
  [string]$SkillsDir   = (Join-Path $env:USERPROFILE ".claude\skills"),
  [string]$ProjetosDir = (Join-Path $env:USERPROFILE "OneDrive\Data_Science\Projetos"),
  [switch]$HookGlobal
)

$ErrorActionPreference = "Stop"
$src = $PSScriptRoot   # ...\<repo>\.claude

Write-Host "===================================================================" -ForegroundColor Cyan
Write-Host " Instalando skills pessoais (Claude Code)" -ForegroundColor Cyan
Write-Host "===================================================================" -ForegroundColor Cyan

# ---- 1) copiar as skills -> %USERPROFILE%\.claude\skills ----
function Copy-Skill([string]$nome) {
  $de = Join-Path $src ("skills\" + $nome)
  if (-not (Test-Path -LiteralPath $de)) {
    Write-Host "   [!] skill origem nao encontrada: $de" -ForegroundColor Red
    return
  }
  New-Item -ItemType Directory -Force -Path $SkillsDir | Out-Null
  $para = Join-Path $SkillsDir $nome
  if (Test-Path -LiteralPath $para) { Remove-Item -LiteralPath $para -Recurse -Force }
  Copy-Item -LiteralPath $de -Destination $para -Recurse -Force
  Write-Host "   [skill] $nome  ->  $para" -ForegroundColor Green
}

Copy-Skill "salvar_contexto"
Copy-Skill "salvar_contexto_compactar"

# ---- 2) copiar o manual PDF -> Projetos ----
$pdf = Join-Path $src "manual_skills.pdf"
if (Test-Path -LiteralPath $pdf) {
  if (Test-Path -LiteralPath $ProjetosDir) {
    Copy-Item -LiteralPath $pdf -Destination (Join-Path $ProjetosDir "manual_skills.pdf") -Force
    Write-Host "   [pdf]   manual_skills.pdf  ->  $ProjetosDir" -ForegroundColor Green
  } else {
    Write-Host "   [!] pasta de destino do PDF nao existe: $ProjetosDir" -ForegroundColor Yellow
  }
} else {
  Write-Host "   [!] manual_skills.pdf nao encontrado ao lado do script." -ForegroundColor Yellow
}

# ---- 3) (opcional) hook global de auto-load ----
if ($HookGlobal) {
  Write-Host " Instalando hook global de auto-load..." -ForegroundColor Cyan
  $hooksDir = Join-Path $env:USERPROFILE ".claude\hooks"
  New-Item -ItemType Directory -Force -Path $hooksDir | Out-Null
  $helperSrc = Join-Path $src "hooks\carregar_contexto.py"
  $helperDst = Join-Path $hooksDir "carregar_contexto.py"
  Copy-Item -LiteralPath $helperSrc -Destination $helperDst -Force
  # Resolve o caminho absoluto do python agora, para o hook nao depender do PATH
  # do processo que o Claude Code usar ao disparar o SessionStart.
  $py = (Get-Command python -ErrorAction SilentlyContinue).Source
  if (-not $py) { $py = (Get-Command python3 -ErrorAction SilentlyContinue).Source }
  if (-not $py) { $py = "python" }
  $cmd = '"' + $py + '" "' + $helperDst + '"'

  $settingsPath = Join-Path $env:USERPROFILE ".claude\settings.json"
  $entry = [pscustomobject]@{
    matcher = ""
    hooks   = @([pscustomobject]@{ type = "command"; command = $cmd })
  }

  if (-not (Test-Path -LiteralPath $settingsPath)) {
    $obj = [pscustomobject]@{ hooks = [pscustomobject]@{ SessionStart = @($entry) } }
    $obj | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $settingsPath -Encoding UTF8
    Write-Host "   [hook]  settings.json global criado com SessionStart." -ForegroundColor Green
  } else {
    $raw = Get-Content -LiteralPath $settingsPath -Raw
    if ($raw -like "*carregar_contexto.py*") {
      Write-Host "   [hook]  ja presente no settings.json global (idempotente)." -ForegroundColor DarkGray
    } else {
      Copy-Item -LiteralPath $settingsPath -Destination ($settingsPath + ".bak") -Force
      try {
        $json = $raw | ConvertFrom-Json
        if (-not ($json.PSObject.Properties.Name -contains "hooks")) {
          $json | Add-Member -NotePropertyName hooks -NotePropertyValue ([pscustomobject]@{}) -Force
        }
        if (-not ($json.hooks.PSObject.Properties.Name -contains "SessionStart")) {
          $json.hooks | Add-Member -NotePropertyName SessionStart -NotePropertyValue @() -Force
        }
        $json.hooks.SessionStart = @(@($json.hooks.SessionStart) + $entry)
        $json | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $settingsPath -Encoding UTF8
        Write-Host "   [hook]  SessionStart adicionado ao settings.json global (backup .bak)." -ForegroundColor Green
      } catch {
        Write-Host "   [!] Nao consegui mesclar o settings.json global automaticamente." -ForegroundColor Red
        Write-Host "       Backup salvo em: $settingsPath.bak" -ForegroundColor Yellow
        Write-Host "       Adicione manualmente este hook em SessionStart:" -ForegroundColor Yellow
        Write-Host "         command: $cmd" -ForegroundColor Yellow
      }
    }
  }
  Write-Host "   Hook global ativo: vale em TODOS os projetos" -ForegroundColor DarkGray
  Write-Host "   (le docs/historico/CONTEXTO_SESSAO.md relativo a pasta de cada projeto)." -ForegroundColor DarkGray
  Write-Host "   Comando do hook: $cmd" -ForegroundColor DarkGray
}

Write-Host ""
Write-Host " Concluido. Abra/renicie o Claude Code para as skills aparecerem." -ForegroundColor Green
Write-Host " Chame com:  /salvar_contexto   ou   /salvar_contexto_compactar" -ForegroundColor Green
