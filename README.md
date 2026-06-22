# NQ-Docking Suite v2 — Pipeline Computacional Antiparasitário

**Projeto:** Maciel et al. — *T. cruzi* + *Leishmania amazonensis*
**Composto base:** 2-Tiocianato-1,4-Naftoquinona (2-SCN-NQ)
**SMILES:** `O=C1C=C(SC#N)C(=O)c2ccccc21`
**Fórmula:** C₁₁H₅NO₂S | **MW:** 215.23 g/mol
**Manuscrito-alvo:** ACS Infectious Diseases (Q1, IF 3.8)

---

## Início Rápido (3 comandos)

```bat
:: 1. Ativar env datascience (conda)
conda activate datascience

:: 2. Rodar o dashboard (porta 5001)
python dashboard_v2.py

:: 3. MD real com GPU (seção MD do dashboard) — requer admin para power-cap
.\run_md_safe.ps1 -Ns 20 -Tag TcTR_inhibitor -PowerLimit 120
```

Abrir no navegador: **http://localhost:5001**

> ⚡ **GPU power-cap obrigatório:** o MD a 100% de GPU (RTX 3060, 161 W) pode
> estourar o no-break (Kernel-Power Event 41). Sempre rodar o dashboard como
> **Administrador** para o cap funcionar, ou aplicar manualmente:
> `nvidia-smi -pl 120` (PowerShell admin).

---

## Módulos e Arquivos Principais

| Arquivo | Função |
|---|---|
| `dashboard_v2.py` | Interface web — todas as abas (:5001) |
| `md_simulation.py` | MD real OpenMM (CUDA/OpenCL) + resume pós-queda |
| `docking_v2.py` | Pipeline AutoDock Vina (prep + QC + docking) |
| `compound_db.py` | Banco de moléculas; **todas as saídas em `figures/`** |
| `dual_design_pipeline.py` | Drug design dual (subversivo + inibidor) — RDKit local |
| `admet_module.py` | ADMET: Lipinski, Caco-2, QED, BOILED-Egg |
| `cyclodex_module.py` | Complexação β-ciclodextrina |
| `analysis.py` | Figuras docking: heatmap, seletividade |
| `vina.exe` | AutoDock Vina 1.2.7 |
| `docs/PROTOCOLO_CIENTIFICO.md` | Parâmetros científicos, alvos, narrativa mecanística |

---

## Estrutura de Diretórios (simplificada — 2026-05-29)

```
docking/
├── dashboard_v2.py
├── md_simulation.py
├── docking_v2.py
├── compound_db.py
├── dual_design_pipeline.py
├── vina.exe
├── run_md_safe.ps1            ← lançador MD seguro (power-cap)
├── subprojeto/                ← padrão guarda-chuva (registry + SUB01 + pipeline SVG)
│   ├── subprojeto_config.json
│   └── SUB01/  (RESUMO_SUB01.md + status.json + pipeline_sub01.svg)
│
└── docking_project/
    ├── proteins/              ← receptores PDB (compartilhados)
    ├── ligand/2-SCN-NQ/       ← ligante preparado
    ├── figures/2-SCN-NQ/      ← TODAS as saídas por sessão (timestamp)
    │   ├── <sess>/            ← figures + docking + logs por sessão
    │   └── md_clean/          ← trajetórias MD (DCD, CSV, CHK, .mode)
    │       ├── md_TcTR_inhibitor.*   ← Fase A (modo inibidor)
    │       ├── md_TcTR_subversive.*  ← Fase B (substrato subversivo FAD)
    │       └── smoke_test_fad/       ← XMLs FAD/LIG, poses iniciais
    └── results/2-SCN-NQ/      ← results.json acumulativo
```

> **Nota:** A pasta `docking/<compound>/` foi **removida** em 2026-05-29 — era vazia.
> Todas as saídas (resultados de docking, figuras, logs) estão em `figures/<compound>/`.

---

## Fluxo Recomendado (2026)

1. **Docking** (aba ⚙ Pipeline) → alvos TcTR/LiTR/hGR + painel completo → Vina 1.2.7
2. **ADMET** (aba 💊) → SMILES → propriedades + BOILED-Egg
3. **MD** (aba ⚡ MD Simulation):
   - Modo **Inibidor** (`tag TcTR_inhibitor`) → confirma ligação na fenda de entrada (~9 Å do FAD)
   - Modo **Subversivo** (`tag TcTR_subversive`) → testa pose na si-face do FAD (3,60 Å da isoaloxazina)
   - Comparação A×B via `analyze_binding_mode.py` (PBC-correto, mdtraj)
4. **Drug Design Dual** (aba 🧬 Drug Design → Dual) → 5 estruturas maiores (subversivo+inibidor) via RDKit local
5. **Toxicidade** (aba ☢) → auto-carrega IC50/CC50/ED50 do 2-SCN-NQ → IS/IT calculados
6. **Paper AI** (aba ✎) → texto para o manuscrito

---

## MD — Dois Modos de Hipótese Mecanística

| Modo | Tag sugerido | Pose inicial | Significado |
|---|---|---|---|
| **Inibidor clássico** | `TcTR_inhibitor` | `lig_TcTR.pdb` (periférica) | Liga na fenda de entrada; bloqueia acesso da tripanotiona |
| **Substrato subversivo** | `TcTR_subversive` | `lig_TcTR_fadsite.pdb` (si-face FAD) | Compite por elétrons do FAD; gera radicais livres |

> **Achado atual (2026-06-20):** o desfecho central é a **pose subversiva (face-_si_ do FAD)**,
> estável em **MD n=3 de ~50 ns** nas duas redutases parasitárias (lig→FAD TcTR 2,26±0,07; LiTR
> 2,43±0,23 Å) e **perdida na hGR humana** (9,28±1,44 Å) — base estrutural da seletividade
> parasito-vs-hospedeiro. A pose **inibidora** (fenda periférica, ~13,7 Å do FAD; controle de 50 ns)
> confirma o contraste. Estatística: ANOVA de Welch p=0,012. Resultados propagados ao manuscrito
> (ACS Infect. Dis.).

---

## Resume pós-queda de energia

O `md_simulation.py` detecta automaticamente o checkpoint (`.chk`) e retoma de onde parou:

```bat
:: Retoma automaticamente (default)
python md_simulation.py --fresh --ns 20 --tag TcTR_inhibitor --mode inhibitor

:: Forçar recomeço do zero (ignora .chk)
python md_simulation.py --fresh --ns 20 --tag TcTR_novo --mode inhibitor --no-resume
```

Cada run grava `md_<tag>.mode` para evitar que o resume de uma run **inibidor** carregue
o checkpoint numa run **subversivo** (e vice-versa). Tags diferentes = runs independentes.

---

## Drug Design Dual (local, sem tokens)

```bat
cd docking/
python dual_design_pipeline.py          # RDKit + Vina (CPU) — gera 5 estruturas + ranking
python dual_design_pipeline.py --no-dock  # apenas RDKit (sem docking, <10s)
python dual_design_pipeline.py --ollama   # + narrativa LLM local
```

Saídas em `docking_project/figures/dual_design/`:
- `dual_design_results.json` — dados estruturados
- `dual_design_report.html` — relatório com ΔG dos dois sítios e imagens 2D

---

## Instalação

### Env datascience (conda — recomendado)
```bat
conda activate datascience
:: já inclui: rdkit, openbabel, openmm 8.2, mdtraj 1.11, flask, numpy, matplotlib
```

### Pip (venv legado)
```bat
pip install flask numpy matplotlib scipy requests rdkit openbabel-wheel biopython mdtraj
```

### AutoDock Vina 1.2.7
Baixar de https://github.com/ccsb-scripps/AutoDock-Vina/releases → renomear para `vina.exe` → mesma pasta dos scripts.

---

## Solução de Problemas

| Problema | Causa | Solução |
|---|---|---|
| Power-cap falha no dashboard | Dashboard sem admin | Rodar como Administrador **ou** `nvidia-smi -pl 120` (admin PS) |
| MD cai na equilibração | GPU 161 W sem cap, no-break estoura | Power-cap antes do warm-up |
| `checkpoint wrong particles` | Tag igual, modo diferente | Usar tag distinto p/ cada modo; a trava de modo agora bloqueia isso |
| `No template for LIG, bonds differ` | CONECT ausente no `_fadsite.pdb` | Regenerar: `python build_fad_siface_pose.py` |
| `pybel: smi not recognised` | BABEL_LIBDIR não setado (conda) | Usar `obabel.exe` CLI (pipeline usa automaticamente) |
| Canvas MD congelado/truncado | `clientWidth=0` quando aba oculta | Corrigido: usa `getBoundingClientRect + parentNode.offsetWidth` |
| `E:\data\docking\docking\` não existe | Pasta removida — era vazia | Saídas agora em `figures/<compound>/<session>/` |
| `No module 'openbabel'` | Não instalado | `conda install -c conda-forge openbabel` |
| Vina não encontrado | Executável ausente | Baixar e colocar na pasta do projeto |
| `0 results generated` | Múltiplas causas | `python diagnostic.py` |

---

## Citações dos Softwares (para seção Métodos)

```
AutoDock Vina v1.2.7:
Eberhardt et al. (2021) J. Chem. Inf. Model. 61:3891–3898.

Open Babel / obabel:
O'Boyle et al. (2011) J. Cheminform. 3:33.

OpenMM v8.2:
Eastman et al. (2017) PLOS Comput. Biol. 13:e1005659.

MDTraj v1.11:
McGibbon et al. (2015) Biophys. J. 109:1528–1532.

RDKit:
Landrum (2023) RDKit: Open-source cheminformatics. https://www.rdkit.org

Protein Data Bank:
Berman et al. (2000) Nucleic Acids Res. 28:235–242.

PyMOL (portátil, renderização 3D):
Schrödinger LLC, The PyMOL Molecular Graphics System, Version 2.5.
```

---

## Hardware

| Componente | Especificação |
|---|---|
| CPU | AMD Ryzen 5 5500 (6c/12t) |
| GPU | NVIDIA GeForce RTX 3060 12 GB |
| RAM | 32 GB DDR4 |
| SO | Windows 11 |
| Python/env | conda datascience (Python 3.11) |
| OpenMM | 8.2.0 (CUDA + OpenCL) |
| Vina | 1.2.7 |

**Tempos MD estimados (sistema ~185k átomos TcTR+FAD+LIG+água):**
- CUDA RTX 3060 (power-cap 120 W): ~14–15 h / 20 ns
- CUDA sem cap (161 W): ~11 h / 20 ns (risco no-break)
- CPU apenas: inviável (>200 h)
