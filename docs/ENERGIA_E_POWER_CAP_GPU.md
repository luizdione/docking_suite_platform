# Energia / Power-cap da GPU para MD — registro e reversão

**Data:** 2026-05-25
**Contexto:** rodadas de MD (`md_clean.py`, CUDA na RTX 3060) provocando
**desligamentos abruptos** do PC e **apito do no-break**.

## Diagnóstico (por que isto foi feito)

Não é bug de software — é **energia/hardware**:

- **Visualizador de Eventos:** `Kernel-Power ID 41, Reason=0` ("reiniciado sem
  desligamento correto / ficou sem energia") em **06/05, 15-17/05 e 24/05** —
  ou seja, **recorrente** e anterior à maior parte do projeto.
- **WHEA-Logger vazio** → sem erro de máquina de CPU/RAM; e **sem BugCheck/tela
  azul** → não é crash de driver.
- **No-break apitando durante a MD** → a carga do PC sob MD (pico GPU+CPU
  ~350-400 W) **excede a capacidade do no-break** (ou a bateria está degradada),
  disparando a proteção → corte de energia → Event 41.
- A MD em si está **saudável** (minimização + warm-up até 300 K + NVT estáveis,
  0 NaN). O que falta é **energia elétrica estável**.

## ⚠️ Carga CONCORRENTE (RNA-seq + MD) — provavel gatilho principal

Há uma **análise de RNA-seq rodando em paralelo** (CPU + RAM intensos). Rodar
RNA-seq (CPU ~100% nos 6c/12t do Ryzen) **junto** com a MD (GPU ~100% na RTX
3060) **soma os dois maiores consumos** → pico combinado de energia bem acima da
MD sozinha → estoura o no-break (apito) e dispara o corte.

**Mitigação nº 1 (de graça, mais eficaz): NÃO rodar os dois ao mesmo tempo.**
Serialize — termine o RNA-seq, depois rode a MD (ou vice-versa). Só isso já
remove o pico combinado e provavelmente elimina os cortes, sem capar nada.

Se PRECISAR sobrepor:
- cap da GPU (abaixo) **+** limitar threads do RNA-seq (deixe 2-4 núcleos livres;
  ex.: `--threads`/`-p` do alinhador/montador em vez de usar todos);
- vigiar a **RAM (32 GB)**: RNA-seq (alinhamento/montagem) pode usar 10-30 GB; se
  RNA-seq + MD passarem de 32 GB → swap/oom (trava por memória, não por energia).

## O que foi feito (paliativo)

Limitar a potência da GPU para reduzir o pico de consumo e tirar carga do
no-break/fonte:

```powershell
# PowerShell como ADMINISTRADOR, antes de rodar a MD:
nvidia-smi -pl 120        # RTX 3060: TGP padrao ~170 W -> limita a 120 W
# (pode tentar 100 se ainda apitar/desligar)
```

Conferir os limites suportados pela placa:
```powershell
nvidia-smi -q -d POWER    # mostra Min/Max/Default Power Limit
```

## Impacto em desempenho (responde: "inviabiliza 50 ns?")

**Não.** A vazão de MD cai pouco com o cap, porque potência escala de forma
sub-linear com clock/desempenho:

| Power limit | Energia | Velocidade aprox. | Risco de corte |
|---|---|---|---|
| 170 W (padrão) | 100% | 100% | **alto** (desliga) |
| 120 W | ~70% | ~80-90% | baixo |
| 100 W | ~59% | ~70-75% | muito baixo |

- 50 ns continua sendo uma rodada **overnight** em qualquer caso; o cap adiciona
  só ~10-25% ao tempo.
- **Uma rodada capada que termina >> uma rodada a 170 W que é cortada.**
- Meça a velocidade real (ns/dia) na coluna `Speed (ns/day)` do CSV de produção
  (`docking_project/figures/md_clean/md_<tag>.csv`) para planejar o tempo.
- Resiliência: produção grava **checkpoint a cada 50k passos**
  (`md_<tag>.chk`) → dá para retomar se cair.

## Como REVERTER

O `nvidia-smi -pl` **não é persistente** — reseta sozinho ao **reiniciar o PC**.
Para reverter manualmente sem reboot:

```powershell
nvidia-smi -pl 170        # volta ao padrao da RTX 3060 (confirme o Default em -q -d POWER)
```

(Se algum dia quiser tornar o cap permanente, seria preciso uma tarefa no
Agendador rodando no boot — NÃO foi feito aqui; o padrão volta ao reiniciar.)

## Correção DURÁVEL (recomendado para produção longa/rotineira)

O cap é paliativo. Para 50 ns e fármacos novos sem sustos, atacar a fonte de energia:

1. **No-break:** trocar a bateria (se velha) ou usar um com **VA/W suficiente**
   para PC com GPU dedicada; OU, em último caso, ligar o PC **direto numa boa
   tomada** (um no-break falhando pode atrapalhar mais que ajudar).
2. **Fonte (PSU):** verificar potência/idade; uma PSU de qualidade com folga
   (≥650 W) elimina o trip sob carga.
3. **Nuvem para a produção:** rodar a fase pesada (produção de ns) em GPU
   alugada (Colab/vast.ai/etc.) — o sistema já está montado e validado; leva-se
   só os arquivos. Mantém este PC para preparo/análise.

## Comando de referência (rodada segura neste PC)

```powershell
nvidia-smi -pl 120
$env:MD_CLEAN_NO_RESTRAINT="1"
python md_clean.py --fresh --ns 0.2 --tag prod_v3 --platform CUDA `
  --fad-xml "docking_project\figures\smoke_test_fad\cofactor_FAD_TcTR_v3.xml"
# Se fechar sem desligar -> escalar --ns gradualmente (1, 5, ... 50),
# monitorando temp/Power no HWiNFO e a coluna Speed no CSV.
```
