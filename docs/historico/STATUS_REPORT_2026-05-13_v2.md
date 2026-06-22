# Status Report — 2026-05-13 v2 (sessão automatizada / validação do fix v2)

**Tarefa:** retomar do limite de créditos. Verificar o estado do `md_clean.py`
após a iteração v2 (parent-anchored H placement, 00:14 UTC) e confirmar
que o pipeline FRESH chega à minimização sem o crash de template do FAD.

---

## 1. O que mudou desde o report anterior (22:33 → 00:25)

Entre o STATUS_REPORT_2026-05-13.md original (22:33) e esta sessão:

| Hora (UTC) | Evento | Arquivo |
|---|---|---|
| 22:33 | Report v1 escrito (Kabsch global, 709 linhas) | `STATUS_REPORT_2026-05-13.md` |
| 23:50 | Run `--tag dryrun_fix` chegou a minimizar | `minimized_dryrun_fix.pdb` |
| ~00:10 | Detectado: H propagado via rotação global cai a 5–10 Å do parent (FAD não é rígido — template e xtal têm RMSD heavy 6–9 Å) | — |
| 00:14 | `md_clean.py` reescrito para v2 (709→749 linhas, +40) | `md_clean.py` |
| 00:22 | Run `--tag dryrun_v2` reconstruiu o sistema | `complex_fresh.pdb` (15 MB) |
| 00:25 | Minimização v2 concluiu | `minimized_dryrun_v2.pdb` (15 MB) |

**Mudança chave (v1 → v2)** em `build_canonical_fad_topology`: cada H
passa a ser ancorado **localmente** ao seu *parent heavy* via a topologia,
em vez de propagado pela rotação global do Kabsch. Comprimento de
ligação X-H é exatamente preservado, e o H fica sempre próximo do átomo
pesado correto na conformação cristalográfica dobrada.

---

## 2. Estado verificado do `md_clean.py` (sessão automática 2026-05-13)

```
AST                         OK
Linhas                      749
Funções top-level (9)       build_canonical_fad_topology, build_system_fast,
                            build_system_fresh, pick_platform,
                            add_position_restraints, diagnose_system,
                            run_md, analyze, main
__main__ guard              OK
Vec3 try/except fallback    presente (linhas 168–171)
addHydrogens(ff, [protein]) ANTES de FAD/LIG (linha 284)
Helper Kabsch + parent-     presente (linhas 132–164)
  anchored H placement
```

Marker "FIX 2026-05-13 v2" confirmado na linha 121:
> *Kabsch global aplica rotação rígida única, MAS o FAD não é rígido…*
> *Solução: para cada H, ancorar localmente no parent heavy.*

---

## 3. Validação programática pura-numpy

Reproduzi a lógica do helper sem OpenMM, lendo os PDBs e CONECT records
direto do disco:

| Caso | Kabsch RMSD heavy | H placement | X-H bond length (Å) |
|---|---|---|---|
| Identidade (template ↔ template-heavy) | 2.6e-16 nm | 33/33 local | erro = 4.4e-16 Å |
| FAD#1 real (`cofactor_FAD_0_xtal.pdb`) | 8.71 Å | 33/33 local, 0 fallback | mean 1.064, range 0.974–1.105 |
| FAD#2 real (`cofactor_FAD_1_xtal.pdb`) | 6.63 Å | 33/33 local, 0 fallback | mean 1.064, range 0.974–1.105 |

O Kabsch RMSD alto (6–9 Å) é **esperado**: o template `cofactor_FAD_TcTR.pdb`
está numa conformação genérica estendida, enquanto os xtal estão nas
conformações dobradas do sítio ativo de TcTR. O importante é que os
comprimentos X-H caem todos na faixa fisiológica:

- C-H esperado ~1.09 Å → observado: vários valores em 1.07–1.10 Å ✓
- N-H esperado ~1.01 Å → observado: dois valores ~1.00 Å ✓
- O-H esperado ~0.97 Å → observado: vários valores ~0.97–0.98 Å ✓

→ A geometria local de grupos rígidos (CH3, CH2, NH2, OH) é preservada
exatamente em todos os 33 hidrogênios.

---

## 4. Inventário do sistema minimizado (v2)

`docking_project/figures/md_clean/minimized_dryrun_v2.pdb`, escrito
2026-05-14 00:25, 15.1 MB, **184.519 átomos**:

| Componente | Quantidade | Confere? |
|---|---|---|
| Proteína (TcTR dímero) | 14.942 átomos (969 resíduos) | ✓ |
| FAD (×2) | 172 átomos (86 canônicos cada) | ✓ |
| LIG (2-SCN-NQ) | 36 átomos | ✓ |
| HOH | 169.056 átomos (~56.352 moléculas) | ✓ |
| Box (Å) | 128.9 × 129.0 × 127.7 | ✓ (padding 1 nm) |

Nenhum NaN no PDB; máx. dist. do centro = 106 Å (compatível com box
~129 Å sem blow-up). **O createSystem + LocalEnergyMinimizer passaram
sem ValueError.**

Observação: ambos os FADs aparecem com `resnum=1` no PDB pós-Modeller,
mas em **chains diferentes** (o Modeller renumera por chain). Para
analyze() que separa FAD pelo `resname=FAD`, isso é OK; análises
per-FAD precisam slice por chain ou range de atom index.

---

## 5. O que ainda não rodou

Nenhum desses arquivos existe ainda:

- `md_dryrun_v2.dcd` (trajetória de produção)
- `md_dryrun_v2_rmsd_rmsf.png` (análise pós)
- `md_dryrun_v2_traj.dcd` em `docking_project/figures/`

→ O run de 00:22–00:25 produziu o sistema + minimização, **mas a
produção MD não foi disparada** (ou foi interrompida antes do warm-up,
provavelmente pelo limite de créditos da sessão anterior).

---

## 6. Próximo passo recomendado para o usuário

O bug original (`No template found for residue 970 (FAD)`) está
**resolvido**. O sistema reconstrói, hidrogena, solvata e minimiza
end-to-end com os FADs canônicos. Falta só rodar a produção.

```powershell
# Limpar caches + processos travados
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
Get-ChildItem -Path 'C:\Users\luizd\OneDrive\Data_Science\Projetos\docking' `
  -Recurse -Include 'md_clean*.pyc','__pycache__' -ErrorAction SilentlyContinue |
  Remove-Item -Force -Recurse

conda activate datascience
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking

# Sanidade AST
python -c "import ast; ast.parse(open('md_clean.py',encoding='utf-8').read()); print('md_clean OK')"

# Rodar produção curta (1 ns ≈ 3-5 min em RTX 3060)
python md_clean.py --fresh --ns 1 --tag dryrun_v2_prod --platform CUDA
```

Logs esperados na produção (após reconstrução idêntica à de 00:25):

```
[OpenMM] warm-up T=300K dt=2fs HBonds...
[OpenMM] Producao OK em ~3 min
RMSD/RMSF salvo em md_dryrun_v2_prod_rmsd_rmsf.png
DONE.
```

Se a produção rodar, abrir o caminho para a sequência de 5 ns / 50 ns
de produção do screening de 2-SCN-NQ vs TcTR (Trypanosoma cruzi
trypanothione reductase) que o projeto vinha perseguindo desde
abril/2026.

---

## 7. Cenários de falha residual a observar

| Sintoma | Causa possível | Mitigação |
|---|---|---|
| `Particle coordinate is NaN` no warm-up | H ancorado mas overlap com proteína vizinha após Modeller mesclar topologias | `LocalEnergyMinimizer` já roda (default em `run_md`); minimização v2 não NaN-ou; se ocorrer só em produção, reintroduzir warm-up dt=0.5fs |
| Run trava no `analyze()` mas DCD escrito | Pós-processamento (mdtraj/MDAnalysis) com versão incompatível | Inspecionar `md_dryrun_v2_prod_rmsd_rmsf.png`; rodar `python md_clean.py --analyze-only --tag dryrun_v2_prod` para refazer só a análise |
| FAD se desloca >5 Å da posição cristalográfica em <1 ns | Sem restraint nos cofatores (default) | Adicionar `add_position_restraints` nos átomos heavy do FAD durante NVT eq (já existe a função, basta wirar no `run_md`) |

---

## 8. Arquivos verificados nesta sessão (não modificados)

- [md_clean.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_clean.py) — 749 linhas, AST OK, helper v2 confirmado, addHydrogens-protein-only confirmado.
- [STATUS_REPORT_2026-05-13.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-13.md) — report v1 (Kabsch global) preservado para histórico.

## 9. Arquivos criados nesta sessão

- [STATUS_REPORT_2026-05-13_v2.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-13_v2.md) — este relatório.

---

## 10. Resumo executivo

**Bug do crash `No template found for residue (FAD)` está resolvido.**
A iteração v2 (00:14) refinou o helper Kabsch para ancorar cada H
localmente no parent heavy via a topologia do template, em vez de
propagar todos via rotação global. Validação pura-numpy mostra que
todos os 33 hidrogênios são colocados via parent-anchor (0 fallback)
e comprimentos X-H ficam em 0.97–1.10 Å (C-H/N-H/O-H corretos).

**O run de 00:22–00:25 confirmou end-to-end** que o pipeline FRESH
agora atravessa: PDBFixer → addHydrogens(protein) → 2× FAD canônicos
via Kabsch → LIG → addSolvent(TIP3P+0.15M NaCl) → createSystem(PME) →
LocalEnergyMinimizer, produzindo um sistema de **184.519 átomos** sem
NaN e com inventário correto (TcTR dímero + 2 FAD + 1 LIG + 56k HOH).

**Bloqueio remanescente: produção MD não chegou a disparar** (ou foi
interrompida pelo limite de créditos). Comando para o usuário:
`python md_clean.py --fresh --ns 1 --tag dryrun_v2_prod --platform CUDA`
— ~3 min em RTX 3060 — deve produzir o primeiro DCD + PNG válido
do pipeline limpo.
