<#
  run_md_safe.ps1  -  Rodada de MD "segura" (TcTR + FAD + 2-SCN-NQ)
  Preparado em 2026-05-26 por automacao.  *** NAO rode com RNA-seq ativo. ***

  Por que: rodar RNA-seq (CPU 100%) + MD (GPU 100%) ao mesmo tempo soma os picos
  de consumo e estoura o no-break -> Kernel-Power Event 41 -> o PC desliga.
  Ver docs/ENERGIA_E_POWER_CAP_GPU.md.

  Sequencia automatizada:
    [0] Confirma que o RNA-seq terminou (serializar e a mitigacao no 1)
    [1] Prepara o env 'datascience' (python + DLLs CUDA no PATH, igual ao dashboard)
    [2] PORTAO: validate_ff.py no FAD e no ligante -> ABORTA se reprovar
    [3a] BLINDAGEM CONTRA STANDBY/HIBERNATE (powercfg + SetThreadExecutionState)
    [3b] Power-cap da GPU (nvidia-smi -pl) para nao estourar o no-break
    [4] Roda md_simulation.py (--fresh)
    [5a] Reverte o power-cap para 170 W
    [5b] Restaura timeouts default de sleep/hibernate

  Uso (PowerShell, DE PREFERENCIA COMO ADMINISTRADOR para o -pl funcionar):
    cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
    .\run_md_safe.ps1                    # 20 ns, tag prod_v5, cap 170 W
    .\run_md_safe.ps1 -Ns 5 -Tag teste5  # rodada curta de sanidade
    .\run_md_safe.ps1 -PowerLimit 100    # cap mais agressivo (PC apita menos)
    .\run_md_safe.ps1 -SkipConfirm       # pula a pergunta do RNA-seq
#>
[CmdletBinding()]
param(
  [double]$Ns = 20,
  [string]$Tag = "prod_v5",
  [int]$PowerLimit = 170,
  [ValidateSet("CUDA","OpenCL","CPU","AUTO")][string]$Platform = "CUDA",
  [switch]$SkipConfirm
)

$ErrorActionPreference = "Continue"
$Proj     = "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
$SmokeDir = Join-Path $Proj "docking_project\figures\smoke_test_fad"
$FadXml   = Join-Path $SmokeDir "cofactor_FAD_TcTR.xml"
$FadPdb   = Join-Path $SmokeDir "cofactor_FAD_TcTR.pdb"
$LigXml   = Join-Path $SmokeDir "lig_TcTR_smirnoff.xml"

Write-Host "===================================================================" -ForegroundColor Cyan
Write-Host " Rodada de MD segura | ns=$Ns  tag=$Tag  plataforma=$Platform  cap=$PowerLimit W" -ForegroundColor Cyan
Write-Host "===================================================================" -ForegroundColor Cyan
Set-Location $Proj

# ---------------------------------------------------------------- [0] RNA-seq
if (-not $SkipConfirm) {
  $ans = Read-Host " [0] O RNA-seq JA TERMINOU? MD + RNA-seq juntos podem cortar a energia. Continuar? (S/N)"
  if ($ans -notmatch '^(s|sim|y|yes)$') {
    Write-Host " Abortado: finalize o RNA-seq antes (ou use -SkipConfirm se ja terminou)." -ForegroundColor Yellow
    exit 3
  }
}

# ---------------------------------------------------------------- [1] Env
$EnvRoot = Join-Path $env:USERPROFILE "miniforge3\envs\datascience"
$Py = Join-Path $EnvRoot "python.exe"
if (-not (Test-Path $Py)) {
  Write-Host " [ERRO] python.exe nao encontrado em $EnvRoot" -ForegroundColor Red
  Write-Host "        Ajuste `$EnvRoot no topo do script, ou rode 'conda activate datascience' e troque `$Py por 'python'." -ForegroundColor Red
  exit 2
}
$env:PATH = (@(
  $EnvRoot,
  (Join-Path $EnvRoot "Library\bin"),          # cudart64_*.dll, cufft, etc.
  (Join-Path $EnvRoot "Library\mingw-w64\bin"),
  (Join-Path $EnvRoot "Library\usr\bin"),
  (Join-Path $EnvRoot "Scripts")
) -join ";") + ";" + $env:PATH
$env:CONDA_PREFIX = $EnvRoot
$env:CONDA_DEFAULT_ENV = "datascience"
Write-Host " [1] Env 'datascience' pronto: $Py" -ForegroundColor Green

# ---------------------------------------------------------------- [2] Gate FF
Write-Host " [2] Validando force-fields (gate validate_ff.py)..." -ForegroundColor Cyan
& $Py validate_ff.py $FadXml --pdb $FadPdb
if ($LASTEXITCODE -ne 0) {
  Write-Host " [ABORTADO] FAD reprovado no validate_ff (exit=$LASTEXITCODE). NAO vou rodar a MD." -ForegroundColor Red
  exit 1
}
& $Py validate_ff.py $LigXml
if ($LASTEXITCODE -ne 0) {
  Write-Host " [ABORTADO] Ligante reprovado no validate_ff (exit=$LASTEXITCODE). NAO vou rodar a MD." -ForegroundColor Red
  exit 1
}
Write-Host " [2] OK: FAD e ligante PASS." -ForegroundColor Green

# ---------------------------------------------------------- [3a] Anti-standby
# Por que: standby/hibernate suspendem o processo do OpenMM. Sintoma: a
# minimizacao "leva" 10 h em vez de 4 min (slowdown 150x = wall-clock dormindo).
# Aplicamos DOIS niveis de protecao:
#   (i)  powercfg /change ... timeout-ac 0   -> sistema nao dorme/hiberna
#   (ii) SetThreadExecutionState (Win32)     -> trava ate este PS fechar
# Os defaults sao restaurados em [5b].
Write-Host " [3a] Blindando contra standby/hibernate..." -ForegroundColor Cyan
# (i) salvar defaults para restaurar depois
$prevStandby   = (powercfg /query SCHEME_CURRENT SUB_SLEEP STANDBYIDLE       2>$null | Select-String -Pattern 'Current AC Power Setting Index:\s*0x([0-9a-fA-F]+)' | ForEach-Object { [Convert]::ToInt32($_.Matches[0].Groups[1].Value,16) }) 2>$null
$prevHibernate = (powercfg /query SCHEME_CURRENT SUB_SLEEP HIBERNATEIDLE     2>$null | Select-String -Pattern 'Current AC Power Setting Index:\s*0x([0-9a-fA-F]+)' | ForEach-Object { [Convert]::ToInt32($_.Matches[0].Groups[1].Value,16) }) 2>$null
$prevDisk      = (powercfg /query SCHEME_CURRENT SUB_DISK   DISKIDLE          2>$null | Select-String -Pattern 'Current AC Power Setting Index:\s*0x([0-9a-fA-F]+)' | ForEach-Object { [Convert]::ToInt32($_.Matches[0].Groups[1].Value,16) }) 2>$null
if (-not $prevStandby)   { $prevStandby   = 1800 }   # default Windows = 30 min
if (-not $prevHibernate) { $prevHibernate = 1800 }
if (-not $prevDisk)      { $prevDisk      = 1200 }   # default disk = 20 min
powercfg /change standby-timeout-ac   0 | Out-Null
powercfg /change hibernate-timeout-ac 0 | Out-Null
powercfg /change disk-timeout-ac      0 | Out-Null
# (ii) SetThreadExecutionState: ES_CONTINUOUS|ES_SYSTEM_REQUIRED|ES_AWAYMODE_REQUIRED
$AntiSleep = Add-Type -PassThru -Name AntiSleep_$([guid]::NewGuid().ToString('N')) -Namespace Win32 -MemberDefinition @'
[System.Runtime.InteropServices.DllImport("kernel32.dll")]
public static extern uint SetThreadExecutionState(uint esFlags);
'@
# 0x80000000 ES_CONTINUOUS | 0x00000001 ES_SYSTEM_REQUIRED | 0x00000040 ES_AWAYMODE_REQUIRED
$null = $AntiSleep::SetThreadExecutionState(0x80000041)
Write-Host "      OK: standby/hibernate/disk-timeout = 0 (AC); SetThreadExecutionState=ES_SYSTEM_REQUIRED|ES_AWAYMODE_REQUIRED" -ForegroundColor Green
Write-Host "      Defaults salvos: standby=$prevStandby s, hibernate=$prevHibernate s, disk=$prevDisk s" -ForegroundColor DarkGray

# ---------------------------------------------------------- [3b] Power-cap GPU
Write-Host " [3b] Aplicando power-cap de $PowerLimit W na GPU..." -ForegroundColor Cyan
$capOk = $false
try {
  nvidia-smi -pl $PowerLimit
  if ($LASTEXITCODE -eq 0) { $capOk = $true }
} catch { }
if (-not $capOk) {
  Write-Host " [AVISO] Power-cap nao aplicado (rode o PowerShell COMO ADMINISTRADOR)." -ForegroundColor Yellow
  Write-Host "         Sem cap, o risco de corte aumenta. Confirme que o RNA-seq esta MESMO parado." -ForegroundColor Yellow
  if (-not $SkipConfirm) {
    $c2 = Read-Host "         Continuar a MD SEM o power-cap? (S/N)"
    if ($c2 -notmatch '^(s|sim|y|yes)$') { Write-Host " Abortado a pedido." -ForegroundColor Yellow; exit 4 }
  }
}

# ---------------------------------------------------------------- [4] MD
Write-Host " [4] md_simulation.py --fresh --ns $Ns --tag $Tag --platform $Platform --allow-fallback" -ForegroundColor Cyan
Write-Host "     Acompanhe ao vivo no dashboard (aba 'MD Simulation' -> graficos de T/energia, tag=$Tag)." -ForegroundColor DarkGray
& $Py md_simulation.py --fresh --ns $Ns --tag $Tag --platform $Platform --allow-fallback
$MdExit = $LASTEXITCODE

# ---------------------------------------------------------- [5a] Reverter cap
if ($capOk) {
  Write-Host " [5a] Revertendo power-cap para 170 W (padrao RTX 3060)..." -ForegroundColor Cyan
  nvidia-smi -pl 170 | Out-Null
}

# ---------------------------------------------------------- [5b] Restaurar sleep
Write-Host " [5b] Restaurando timeouts default de sleep/hibernate..." -ForegroundColor Cyan
powercfg /change standby-timeout-ac   $prevStandby    | Out-Null
powercfg /change hibernate-timeout-ac $prevHibernate  | Out-Null
powercfg /change disk-timeout-ac      $prevDisk       | Out-Null
# Libera o flag de execucao continua
$null = $AntiSleep::SetThreadExecutionState(0x80000000)  # ES_CONTINUOUS sozinho = limpa
Write-Host "      OK: standby=$prevStandby s, hibernate=$prevHibernate s, disk=$prevDisk s" -ForegroundColor Green

if ($MdExit -eq 0) {
  Write-Host " CONCLUIDO (exit=0). Para a analise/figuras: $Py md_simulation.py --analyze-only --tag $Tag" -ForegroundColor Green
} else {
  Write-Host " MD terminou com exit=$MdExit. Veja docking_project\figures\md_clean\md_clean.log" -ForegroundColor Yellow
}
exit $MdExit
