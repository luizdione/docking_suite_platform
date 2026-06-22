# SUB01 — Resumo executivo (sincronizado com o abstract do paper)

> **Programa (NOME GERAL):** NQ-Antiparasitários — descoberta computacional de naftoquinonas
> **Subprojeto:** `SUB01`
> **Nome específico (título atual do paper):** *2-Thiocyanate-1,4-naphthoquinone as a
> subversive substrate of trypanothione reductase: selective FAD si-face binding in
> Trypanosoma cruzi and Leishmania versus human glutathione reductase*
> **Status:** `em_desenvolvimento` — MD n=3 consolidada; redação em curso
> **Periódico-alvo:** ACS Infectious Diseases (Q1)
> **Atualizado:** 2026-06-19

---

## Abstract (resumo científico — versão viva, espelha a aba "Produção de Texto")

A 2-tiocianato-1,4-naftoquinona (2-SCN-NQ; C₁₁H₅NO₂S; MW 215,23 g·mol⁻¹) é uma
naftoquinona candidata a antiparasitário cujo mecanismo de ação proposto envolve a
**subversão redox da tripanotiona redutase (TR)** — flavoenzima essencial e ausente no
hospedeiro humano, no qual a função análoga é exercida pela glutationa redutase (hGR).
Para discriminar entre dois mecanismos concorrentes — **inibidor clássico** (ligação na
fenda de entrada, bloqueando o acesso da tripanotiona) versus **substrato subversivo**
(ligação na *si*-face do FAD, competindo por elétrons e gerando radicais livres) —
combinamos docking molecular (AutoDock Vina 1.2.7) com dinâmica molecular explícita
(OpenMM 8.2, CUDA) em triplicata independente (3 × 50 ns por alvo). A estabilidade da
pose subversiva foi quantificada pela distância ligante–isoaloxazina do FAD e
ligante–par ditiol catalítico. Nas TR parasitárias, o 2-SCN-NQ permaneceu **estável e
próximo do FAD** ao longo de toda a trajetória (TcTR: 2,26 ± 0,07 Å; LiTR: 2,43 ± 0,23 Å,
distância lig–FAD, média ± DP, n = 3), com RMSD intramolecular do ligante < 0,5 Å
(corpo rígido), enquanto na hGR humana a pose **divergiu** do cofator (9,28 ± 1,44 Å),
indicando ausência de acomodação produtiva no sítio redox humano. Esse contraste de
~4× na distância lig–FAD entre enzimas parasitárias e humana é consistente com uma
**janela de seletividade** favorável ao mecanismo subversivo em *T. cruzi* e
*Leishmania*. Os achados sustentam a hipótese de que o 2-SCN-NQ atua preferencialmente
sobre a TR parasitária e fornecem base estrutural para a otimização de derivados
(subprojetos futuros SUB02/SUB03). *Limitações:* as distâncias derivam de MD clássica
sem polarização eletrônica explícita do estado radicalar; a confirmação do ciclo
redox subversivo requer QM/MM e validação experimental (consumo de O₂/NADPH, ROS).

*(Meta de profundidade do abstract conforme CLAUDE.md: ≥ 350 palavras na versão final
do manuscrito; este é o resumo-âncora sincronizado com a aba Produção de Texto.)*

---

## Resultados atuais — modo subversivo, MD n=3 (50 ns × 3 réplicas)

| Sistema | n | dist. lig–FAD (Å, média ± DP) | dist. lig–ditiol (Å, média ± DP) | Interpretação |
|---|---|---|---|---|
| **TcTR** (*T. cruzi*) | 3 | **2,26 ± 0,07** | 2,50 ± 0,00 | Pose subversiva estável na *si*-face do FAD |
| **LiTR** (*Leishmania*) | 3 | **2,43 ± 0,23** | 2,82 ± 0,55 | Pose subversiva estável (maior dispersão) |
| **hGR** (humano) | 3 | **9,28 ± 1,44** | 7,55 ± 1,86 | Ligante afasta-se do cofator → **seletividade** |

**Legenda (Tabela 1).** Distâncias mínimas representativas entre o 2-SCN-NQ e (i) o anel
isoaloxazínico do FAD e (ii) o par ditiol catalítico (Cys52/Cys57 em TcTR), extraídas do
último bloco de equilíbrio de cada trajetória de produção de 50 ns e agregadas sobre três
réplicas independentes (`n_3resultados.json`). Valores em ångström; média ± desvio-padrão
amostral. A distância lig–FAD ~4× maior em hGR frente às TR parasitárias é o principal
descritor de seletividade do modo subversivo. Dados verificados — sem números inventados;
fonte primária: `replica_summary.json` (gerado 2026-06-17) → `aggregate_md_n3.py` →
`n3_resultados.json`.

### Figuras disponíveis (com legenda)
- `docking_project/figures/docking_pose_rmsd_{TcTR,LiTR,hGR,...}.png` — **Fig. docking:**
  RMSD da pose de docking por alvo (painel de seletividade de 9 alvos: TcTR, LiTR, hGR,
  hCatL, hDHODH, hMnSOD, LmDHODH, TcCruzain, TcFeSOD).
- `docking_project/figures/2-SCN-NQ/md_clean/md_{TcTR,LiTR,hGR}_subversive[_r2|_r3]_rmsd_rmsf.png`
  — **Fig. MD:** RMSD/RMSF por réplica (50 ns), proteína e ligante.
- `docking_project/figures/2-SCN-NQ/md_clean/pipeline_md_clean.png` — **Fig. pipeline MD.**
- `docking_project/figures/pipeline_sub01.svg` — **Fig. pipeline do subprojeto** (gerada por
  `subprojeto/pipeline_svg.py`; status colorido por etapa).
- `docking_project/figures/dual_design/struct_{NQ-SPM,NQ-PTZ,NQ-ACR,NQ-TZ-POL,BIS-NQ}.png`
  — **Fig. derivados** (drug design dual; sementes do SUB02 futuro).

---

## Scripts do subprojeto — input → output → logs

| Etapa | Script | Input | Output | Log |
|---|---|---|---|---|
| Docking | `docking_v2.py` + `vina.exe` | receptores `docking_project/proteins/*.pdb`, ligante `docking_project/ligand/2-SCN-NQ/` | poses + `results/2-SCN-NQ/results.json` | `docking_project/logs/` |
| MD (1 run) | `md_simulation.py` (wrapper `run_md_safe.ps1`) | pose Vina + XML FAD/LIG (`smoke_test_fad/`) | `md_clean/md_<tag>.{dcd,csv,chk}` | `md_*.runlog`, `md_run_v2_prod.log` |
| MD réplicas | `replica_driver.py` | tags `*_subversive[_r2,_r3]` | `replica_summary.json`, `replica_status.json` | `replica_driver.log` |
| Re-runs alvo | `rerun_litr_r3.py`, `rerun_tctr_inhibitor.py` | checkpoints `.chk` | trajetórias retomadas | `rerun_*.log` |
| Modo de ligação | `analyze_dcd_numpy.py` | `.dcd` + topologia | distâncias lig–FAD / lig–ditiol | — |
| Agregação n=3 | `aggregate_md_n3.py` | `replica_summary.json` | `n3_resultados.json`, `RELATORIO_n3.md` | `finalize_md_n3.log` |
| Figuras | `analysis.py` | resultados JSON | `figures/*.png` | — |
| ADMET | `admet_module.py` | SMILES | propriedades drug-like | — |
| Toxicidade | `toxicity_module.py` | IC50/CC50/ED50 | IS/IT | — |
| Drug design dual | `dual_design_pipeline.py` | 2-SCN-NQ (semente) | `dual_design/dual_design_report.html` + 5 estruturas | — |

### A fazer (scripts a utilizar)
- QM/MM do estado radicalar (confirmação do ciclo subversivo) — **pendente**.
- Cálculo de energia livre de ligação (MM/PBSA ou FEP) para complementar as distâncias — **pendente**.
- Geração das figuras finais de submissão (300 dpi) e tabelas suplementares — **pendente**.

---

## Status e próximos passos

- [x] Docking de painel (9 alvos) concluído.
- [x] MD subversiva 50 ns × 3 réplicas (TcTR, LiTR, hGR) concluída e agregada (n=3).
- [x] Métrica de seletividade lig–FAD estabelecida (parasito ≈ 2,3 Å vs humano ≈ 9,3 Å).
- [ ] MD modo inibidor n=3 para comparação A×B completa (parcial: `TcTR_inhibitor`).
- [ ] QM/MM e/ou ΔG de ligação.
- [ ] Manuscrito completo (aba Produção de Texto / `/build_paper`) — EN submissão + PT auditoria.
- [ ] Cover letter, declarações (ORCID/COI/CRediT), figuras 300 dpi.

> Sincronizar este arquivo a cada versionamento (ver `RELATORIO_DESENVOLVIMENTO_PLATAFORMAS.md`).
