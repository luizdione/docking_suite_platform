#Requires -Version 5.1
<#
.SYNOPSIS
  Publica a imagem Docker da plataforma no registry com tag = commit atual.
  Fecha a REGRA DE VERSIONAMENTO: versionar o repo -> versionar o pacote online
  (o RunPod consome a imagem do registry). A tag da imagem == o short SHA do commit,
  entao da para rastrear exatamente qual codigo esta rodando na nuvem.

.DESCRIPTION
  Por padrao BUILDA a imagem a partir do codigo atual (garante imagem == commit),
  taguea com o short SHA de HEAD e com :latest, e da push nos dois. Requer um
  `docker login` previo no registry -- este script NAO manuseia credenciais.

.PARAMETER Owner     Namespace no registry (default: luizdione).
.PARAMETER Registry  Host do registry (default: ghcr.io).
.PARAMETER Image     Nome da imagem (default: docking_suite_platform).
.PARAMETER NoBuild   Nao builda; retaguea+push de uma imagem local ja existente (-LocalTag).
.PARAMETER Force     Empurra mesmo com working tree sujo (a tag pode nao bater com a imagem).

.EXAMPLE
  ./push_image.ps1
.EXAMPLE
  ./push_image.ps1 -NoBuild            # so retag+push da imagem local
.EXAMPLE
  ./push_image.ps1 -Owner luizdione -Registry ghcr.io
#>
[CmdletBinding()]
param(
  [string]$Registry = 'ghcr.io',
  [string]$Owner    = 'luizdione',
  [string]$Image    = 'docking_suite_platform',
  [string]$LocalTag = 'test',
  [switch]$NoBuild,
  [switch]$Force
)

$ErrorActionPreference = 'Stop'
# Docker Desktop: client (1.55) mais novo que o engine (1.44) -> forcar versao.
$env:DOCKER_API_VERSION = '1.44'

# Raiz do repo = pasta deste script (roda de qualquer CWD).
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Root

function Assert-LastExit([string]$what) {
  if ($LASTEXITCODE -ne 0) { throw "FALHOU: $what (exit $LASTEXITCODE)" }
}

# 1) commit atual + estado da arvore ------------------------------------------
$sha = (git rev-parse --short HEAD).Trim(); Assert-LastExit 'git rev-parse'
$dirty = git status --porcelain
if ($dirty -and -not $Force) {
  throw "working tree sujo -- a tag :$sha nao refletiria a imagem. Commite/limpe, ou use -Force."
}
if ($dirty) { Write-Warning "working tree sujo; prosseguindo com -Force (tag :$sha pode nao bater com a imagem)." }

$Repo = "$Registry/$Owner/$Image"
Write-Host "==> alvo: ${Repo}  (tags: $sha, latest)" -ForegroundColor Cyan

# 2) build (default) ou imagem local existente --------------------------------
if ($NoBuild) {
  $src = "${Image}:${LocalTag}"
  Write-Host "==> -NoBuild: usando imagem local $src"
  docker image inspect $src *> $null; Assert-LastExit "imagem local $src nao existe (rode SEM -NoBuild p/ buildar)"
} else {
  $src = "${Image}:${sha}"
  Write-Host "==> docker build -t $src ."
  docker build -t $src .; Assert-LastExit 'docker build'
  Write-Host "==> sanidade: vina --version na imagem"
  docker run --rm --entrypoint vina $src --version; Assert-LastExit 'sanidade vina na imagem'
}

# 3) tag para o registry ------------------------------------------------------
docker tag $src "${Repo}:${sha}";   Assert-LastExit 'docker tag sha'
docker tag $src "${Repo}:latest";   Assert-LastExit 'docker tag latest'

# 4) push (requer docker login previo; falha alto se nao autenticado) ---------
Write-Host "==> docker push ${Repo}:${sha}" -ForegroundColor Cyan
docker push "${Repo}:${sha}";   Assert-LastExit 'docker push sha (voce fez docker login?)'
docker push "${Repo}:latest";   Assert-LastExit 'docker push latest'

Write-Host ""
Write-Host "OK -- publicado:" -ForegroundColor Green
Write-Host "  ${Repo}:${sha}"
Write-Host "  ${Repo}:latest"
Write-Host ""
Write-Host "No RunPod: aponte a imagem para  ${Repo}:latest  (ou :${sha} p/ travar a versao)."
Write-Host "Lembrete: pacote GHCR novo nasce PRIVADO -- torne publico ou cadastre credencial no RunPod."
