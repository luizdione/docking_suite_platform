# Relatório de Status — Diagnóstico Pós-MD 50ns
**Data:** 2026-05-07 | **Atualização tras correção crítica**

---

## 1. Achado Crítico — A MD de 50 ns não validou interação proteína-ligante

A simulação de 50 ns documentada no log do dia rodou com o **PDBQT de preparação** (`2-SCN-NQ.pdbqt`, 1641 bytes — apenas a molécula isolada parametrizada para o Vina), **não** com a pose docked. A linha indicativa no log foi:

```
[Ligante] Preparação: C:\...\docking_project\ligand\2-SCN-NQ\2-SCN-NQ.pdbqt
```

Como o ligante não foi posicionado no bolsão da TcTR, os resultados de RMSD (Backbone 2.16 Å, Ligand 0.86 Å) **não confirmam ancoragem no sítio ativo** — apenas indicam dinâmica da apo-proteína mais uma molécula livre.

**Diagnóstico do `check_docking_status.py`:** os 7 alvos retornam `preparation_ONLY` e nenhum `result_<TID>.pdbqt` existe na máquina. **O docking molecular nunca foi executado** nesta pasta.

---

## 2. Correções Aplicadas no md_module.py

| # | Correção | Linha | Status |
|---|---|---|---|
| 1 | `COFACTOR_SMILES` dict (FAD, FMN, NAD, NDP, NAP, HEM) | 772 | ✓ |
| 2 | `_extract_and_parametrize_cofactors()` — parametriza cofatores via SMIRNOFF/RDKit | 792 | ✓ |
| 3 | Chamada do passo "1b" após carregar receptor | 1021 | ✓ |
| 4 | Whitelist preserva cofatores parametrizados | 1431 | ✓ |
| 5 | Aborto explícito quando `result_<TID>.pdbqt` ausente | 2975 | ✓ |
| 6 | AST OK | — | ✓ (3091 linhas, 26 funções, 1 classe) |

**Backup do arquivo quebrado anterior:** `md_module.py.broken_2026-05-07_backup` (141.412 bytes, AST falha)

**Arquivo atual:** `md_module.py` (150.825 bytes, AST OK)

---

## 3. Como o md_module.py se comporta agora

### 3.1 Quando há pose docked (`result_<TID>.pdbqt`)
```
[Cofator] Detectados: FAD
[OK] Cofator FAD parametrizado e registrado (84 átomos)
[OpenMM] Cofatores parametrizados mantidos: FAD
[OpenMM] Removing non-standard residues: GCG
[Ligante] Pose docked: E:\data\docking\docking\2-SCN-NQ\<sessão>\result_TcTR.pdbqt
... (segue normalmente para 50 ns)
```

### 3.2 Quando falta pose docked
```
[Ligante] Preparacao (sem pose docked): C:\...\ligand\2-SCN-NQ\2-SCN-NQ.pdbqt
[ERRO CRITICO] result_TcTR.pdbqt nao encontrado. MD sem pose docked produz...
[ABORT] Pulando MD para TcTR - faltando pose docked do Vina.
```
A função retorna `{"status": "aborted_no_docked_pose", ...}` — sem desperdício de GPU.

---

## 4. Próximos Passos (em ordem)

### FASE 1 — Rodar o docking (URGENTE, primeiro)
```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
.\run_docking_batch.bat
```
Isso vai:
- Baixar PDBs do RCSB (1BZL, 2JK6, 3KKU, 4DVH, 3MHU, 3GRS, 2XU3)
- Preparar receptores PDBQT com cofatores retidos
- Rodar Vina (exhaustiveness=32) para os 7 alvos
- Salvar `docking_project\docking\result_<TID>.pdbqt` para cada
- ~5-10 min por alvo, total ~1 hora

### FASE 2 — Validar com `check_docking_status.py`
```powershell
python check_docking_status.py
```
Esperado: todos os alvos com `[OK] docked_session` ou `[OK] docked_flat`.

### FASE 3 — Teste curto de MD (10 ns) com pose docked + FAD
```powershell
python md_module.py --targets TcTR --ns 10 --mode real
```
**Verificar no log:**
- `[Cofator] Detectados: FAD` ← novo
- `[OK] Cofator FAD parametrizado e registrado` ← novo
- `[OpenMM] Cofatores parametrizados mantidos: FAD` ← novo
- `[OpenMM] Removing non-standard residues:` deve listar **apenas GCG** e similares (nunca FAD)
- `[Ligante] Pose docked: ...result_TcTR.pdbqt` ← novo (não mais "Preparação")

### FASE 4 — Validação visual no PyMOL
```pymol
load docking_project/figures/2-SCN-NQ/<sess>/complex_TcTR_initial.pdb
load_traj docking_project/figures/2-SCN-NQ/<sess>/traj_TcTR.dcd
select lig, resn LIG
select fad, resn FAD
select pocket, byres polymer within 5 of lig
show sticks, lig or fad
zoom lig, 8
```
Confirmar: ligante no bolsão da TcTR, FAD na sua posição cristalográfica, contatos com Cys52/57 / His461 / Glu466.

### FASE 5 — Produção 50-100 ns (apenas após Fase 4 OK)
```powershell
python md_module.py --targets TcTR LiTR TcCruzain hGR --ns 100 --mode real --lang pt
```

---

## 5. Arquivos relacionados

- [md_module.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_module.py) — 3091 linhas, AST OK, com suporte a FAD
- [check_docking_status.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\check_docking_status.py) — diagnóstico de poses docked
- [run_docking_batch.bat](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\run_docking_batch.bat) — pipeline de docking automatizado
- [md_module.py.broken_2026-05-07_backup](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_module.py.broken_2026-05-07_backup) — backup do arquivo truncado

---

## 6. Lições da execução de hoje

1. **Fallback silencioso é perigoso** — o código antes prosseguia com PDBQT de preparação sem alertar; agora aborta explicitamente
2. **Cofatores essenciais não podem ser strippados** — FAD em flavoenzimas (TcTR/LiTR/hGR) é parte do sítio ativo; sem ele a topologia do bolsão muda
3. **OneDrive sync pode mascarar truncação** — a edição via ferramenta vs. estado em disco divergiu; sempre validar AST após mudanças grandes
4. **Pre-condições devem ser verificadas** — rodar 20h de GPU em uma simulação cientificamente inválida é o pior caso; o aborto explícito custa segundos e poupa horas

---

## 7. Resumo Executivo

A MD de 50 ns **não validou** o complexo proteína-ligante porque o ligante nunca foi colocado no bolsão (faltava pose docked). Corrigi `md_module.py` para: (a) suportar FAD via parametrização SMIRNOFF/RDKit automática, (b) abortar explicitamente quando faltar pose docked, e (c) gerei script auxiliar `check_docking_status.py` que confirmou: o docking nunca rodou nesta pasta. **Próximo passo crítico: executar `run_docking_batch.bat` para gerar as poses docked, depois fazer um teste de MD curto (10 ns) com TcTR antes de qualquer produção longa.**
