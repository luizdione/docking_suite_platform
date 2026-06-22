<#
  prevent_sleep_now.ps1  -  Desabilita standby/hibernate AGORA, sem mexer na MD.

  Use quando uma MD ja esta rodando e voce esqueceu de blindar a sessao contra
  sleep. NAO interrompe processo nenhum: so muda os timeouts do esquema de
  energia ativo. Reverte ao final do script (ou manualmente — ver final do output).

  Uso:
    .\prevent_sleep_now.ps1            # bloqueia standby/hibernate/disk; deixa monitor dormir
    .\prevent_sleep_now.ps1 -Restore   # restaura os defaults do Windows (standby 30 min)

  NAO requer Administrador.
#>
[CmdletBinding()]
param(
  [switch]$Restore
)

if ($Restore) {
  Write-Host " Restaurando timeouts default (standby/hibernate 30 min, disk 20 min)..." -ForegroundColor Cyan
  powercfg /change standby-timeout-ac 30
  powercfg /change hibernate-timeout-ac 30
  powercfg /change disk-timeout-ac    20
  powercfg /change monitor-timeout-ac 15
  Write-Host " OK. Esquema de energia voltou aos defaults." -ForegroundColor Green
  exit 0
}

Write-Host "===================================================================" -ForegroundColor Cyan
Write-Host " Blindando a sessao contra standby/hibernate (AC, ligado na tomada)" -ForegroundColor Cyan
Write-Host "===================================================================" -ForegroundColor Cyan

# 0 = nunca dormir / nunca hibernar / disco nunca desliga
powercfg /change standby-timeout-ac 0   # nao entra em standby (S3 / Modern Standby)
powercfg /change hibernate-timeout-ac 0 # nao hiberna
powercfg /change disk-timeout-ac    0   # disco nao para (evita travar I/O do DCD)
# Monitor PODE dormir (so a tela) -- isso nao afeta a MD:
powercfg /change monitor-timeout-ac 15

Write-Host ""
Write-Host " OK. Enquanto o cabo de forca estiver conectado:" -ForegroundColor Green
Write-Host "   * Standby/Sleep:   DESABILITADO" -ForegroundColor Green
Write-Host "   * Hibernate:       DESABILITADO" -ForegroundColor Green
Write-Host "   * Disk timeout:    DESABILITADO (DCD writer fica seguro)" -ForegroundColor Green
Write-Host "   * Monitor:         pode dormir em 15 min (so a tela, nao afeta MD)" -ForegroundColor DarkGray
Write-Host ""
Write-Host " Estado atual do esquema ativo:" -ForegroundColor Cyan
powercfg /query SCHEME_CURRENT SUB_SLEEP STANDBYIDLE 2>$null | Select-String "Power Setting Index|Configuracao de energia" | Select-Object -First 2

Write-Host ""
Write-Host " Quando a MD terminar e voce quiser RESTAURAR os defaults:" -ForegroundColor Yellow
Write-Host "   .\prevent_sleep_now.ps1 -Restore" -ForegroundColor Yellow
