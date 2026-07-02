# =============================================================================
# sync.ps1 - traz as ultimas mudancas do upstream para ESTE repositorio (bonder)
# e publica no seu GitHub (origin = luizdione/bonder).
#
# Contexto: as atualizacoes do assistente sao empurradas para o repositorio
# 'docking_suite' (subpasta bonder/). Este script copia esse conteudo por cima do
# seu repo standalone e faz commit+push no seu origin. Seus dados locais (data/,
# runs/, inputs/ e *.fasta das iscas) NAO sao tocados.
#
# Uso (dentro da pasta do seu repo bonder):
#   powershell -ExecutionPolicy Bypass -File .\sync.ps1
# =============================================================================
param([string]$Branch = "claude/bonder-proteome-platform-k8m6pa")
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

$tmp = Join-Path $env:TEMP "bonder_upstream"
if (Test-Path $tmp) { Remove-Item -Recurse -Force $tmp }

Write-Host "Baixando upstream (docking_suite @ $Branch)..." -ForegroundColor Cyan
git clone --depth 1 --branch $Branch https://github.com/luizdione/docking_suite $tmp
if (-not (Test-Path (Join-Path $tmp "bonder"))) { throw "Subpasta bonder/ nao encontrada no upstream." }

Write-Host "Aplicando arquivos (preserva .git, data, runs, inputs e *.fasta)..." -ForegroundColor Cyan
robocopy "$tmp\bonder" "." /E /XD .git data runs inputs /XF *.fasta | Out-Null
Remove-Item -Recurse -Force $tmp

git add -A
$changes = git status --porcelain
if ($changes) {
    git commit -m "sync: ultimas mudancas do upstream (docking_suite)"
    git push origin main
    Write-Host "OK: sincronizado e publicado em luizdione/bonder (main)." -ForegroundColor Green
} else {
    Write-Host "Nada novo para sincronizar." -ForegroundColor Yellow
}
