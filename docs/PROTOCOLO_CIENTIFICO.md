# Protocolo Científico — Docking Molecular + Perfil Computacional
## 2-Tiocianato-1,4-Naftoquinona (2-SCN-NQ) vs Tripanossomatídeos
### Manuscrito: Santos et al. — ACS Infectious Diseases

**Composto:** 2-SCN-NQ | SMILES: `O=C1C=C(SC#N)C(=O)c2ccccc21`  
**Fórmula:** C₁₁H₅NO₂S | MW: 215.23 g/mol | Carga: 0 | Ligações rotacionáveis: 1

---

## 1. Alvos Moleculares e Justificativas

### 1.1 Trypanosoma cruzi

| ID | PDB | Alvo | Res. | Cofatores mantidos | Justificativa científica |
|---|---|---|---|---|---|
| TcTR | 1BZL | Tripanotiona Redutase | 2.40 Å | FAD | Ausente em mamíferos; conectada ao mecanismo de ROS; alvo clássico naftoquinonas |
| TcCruzain | 3GRS | Cruzaína (cisteíno-protease) | 1.28 Å | nenhum | Essencial para sobrevivência; co-cristal com inibidor covalente |
| TcFeSOD | 1E7W | Fe-Superóxido Dismutase | 1.50 Å | Fe²⁺ | Defesa antioxidante do parasito; diretamente relacionada aos dados de ROS |

### 1.2 Leishmania spp.

| ID | PDB | Alvo | Res. | Cofatores mantidos | Justificativa científica |
|---|---|---|---|---|---|
| LiTR | 2JK6 | Tripanotiona Redutase | 2.95 Å | FAD | Proxy para *L. amazonensis* (95% identidade); único PDB disponível |
| LmDHODH | 3C61 | Diidroorotato Desidrogenase | 2.40 Å | FMN | Via de pirimidinas; alvo emergente em Leishmania |

> **Nota sobre o proxy LiTR:** Não existe estrutura cristalográfica de TR de *L. amazonensis*. O PDB 2JK6 (*L. infantum*) é o proxy padrão na literatura por 95% de identidade de sequência. Referência: Camargo et al., *Sci. Rep.* 2023, 13:6887; Baiocco et al., *Structure* 2009, 17:1547.

### 1.3 Ortólogos Humanos (controles de seletividade)

| ID | PDB | Alvo | Res. | Cofatores | Contraparte parasitária |
|---|---|---|---|---|---|
| hGR | 3DJJ | Glutationa Redutase | 1.54 Å | FAD | TcTR, LiTR ← **par principal** |
| hCatL | 2XU3 | Catepsina L | 1.40 Å | nenhum | TcCruzain |
| hMnSOD | 1LUV | Mn-Superóxido Dismutase | 1.85 Å | Mn | TcFeSOD |
| hDHODH | 4DVH | DHODH humana | 1.60 Å | FMN | LmDHODH |

**Par principal para seletividade:** TcTR (1BZL) vs hGR (3DJJ)
- Enzimas funcionalmente equivalentes (tripanotiona vs glutationa)
- TR tem sítio ativo mais largo para acomodar a ponte espermidina da tripanotiona
- Diretamente conectadas ao mecanismo de estresse oxidativo observado experimentalmente

---

## 2. Parâmetros Técnicos do Docking

### 2.1 Preparação do Ligante
- Geração 3D: Open Babel (pybel) — SMILES → 3D
- Campo de força: MMFF94
- Otimização: 1000 passos gradientes conjugados
- Adição de hidrogênios: sim
- Formato final: PDBQT (gerado automaticamente pelo `compound_db.py`)

### 2.2 Preparação das Proteínas — O que manter e remover

**REMOVIDO:**
- Moléculas de água (HOH, WAT, DOD)
- Íons de tampão (SO₄, PO₄, CL, GOL, EDO, etc.)
- Ligantes co-cristalizados (usados para calcular o centro do grid, depois removidos)

**MANTIDO (crítico!):**

| Molécula | Razão para manter |
|---|---|
| FAD | Forma o sítio ativo da TR; elétrons: NADPH → FAD → Cys52–Cys57 |
| FMN | Cofator essencial da DHODH |
| NADPH (NAP/NDP) | Cofator que forma parte do sítio de ligação |
| Fe²⁺, Mn²⁺, Zn²⁺, Mg²⁺ | Íons catalíticos essenciais nos sítios ativos |

### 2.3 Grid Box — Estratégia de Centralização

| Prioridade | Estratégia | Quando usar |
|---|---|---|
| 1ª | Centro no ligante co-cristalizado | Sempre que houver co-cristal no PDB |
| 2ª | Centro no íon metálico | Fe-SOD, MnSOD |
| 3ª | Centro nos resíduos do sítio ativo | Estruturas apo |

**Tamanhos de grid:**

| Alvo | Tamanho (Å³) | Justificativa |
|---|---|---|
| TcTR, LiTR, hGR | 25 × 25 × 25 | Sítio largo — acomoda tripanotiona |
| TcCruzain, hCatL | 22 × 22 × 22 | Sítio menor — cisteíno-protease |
| TcFeSOD, hMnSOD | 22 × 22 × 22 | Sítio compacto — metaloproteína |
| DHODH | 25 × 25 × 25 | Sítio do ubiquinol — espaçoso |

### 2.4 Parâmetros AutoDock Vina

```
exhaustiveness : 32  (publicação); 8 (teste rápido); 64 (máxima precisão)
num_modes      : 10
energy_range   : 3 kcal/mol
scoring        : vina (padrão)
```

### 2.5 Conversão PDBQT — Atenção Importante

O Open Babel gera PDBQT em formato de **ligante** (com tags ROOT/BRANCH).
O pipeline `docking_v2.py` remove automaticamente essas tags para formato **receptor** válido.
Se preparar manualmente: usar UCSF Chimera → DockPrep → Save as PDBQT.

---

## 3. Dinâmica Molecular (MD)

### 3.1 Dois Modos de Hipótese Mecanística (decisão experimental, 2026-05-29)

O `md_simulation.py` suporta dois modos de hipótese, selecionados por `--mode`:

| Modo | Tag sugerido | Pose inicial do ligante | Hipótese |
|---|---|---|---|
| **`inhibitor`** (inibidor clássico) | `TcTR_inhibitor` | `lig_TcTR.pdb` — fenda de entrada | Liga no bolso hidrofóbico de entrada (~9 Å do FAD); bloqueia acesso da tripanotiona; equivalente ao sítio tipo mepacrina |
| **`subversive`** (substrato subversivo) | `TcTR_subversive` | `lig_TcTR_fadsite.pdb` — si-face do FAD | Empilha 3,60 Å da isoaloxazina; compete por elétrons do FAD; quinona reduzida → semiquinona → cicla com O₂ → superóxido |

**Achado atual (2026-05-29):** Docking (2 grids independentes) + MD 5 ns (Fase A, RMSD 0,5 Å) convergem para o modo **Inibidor** como sítio estável do 2-SCN-NQ. Fase B (subversivo, si-face FAD) em andamento. A comparação A×B definirá a narrativa mecanística final do paper.

> **Importante:** cada modo usa **TAG DISTINTO** para evitar contaminação de checkpoint.
> O script grava `md_<tag>.mode` e rejeita resume de modo errado com erro claro.

### 3.2 Parâmetros recomendados

| Parâmetro | Valor atual (md_simulation.py) |
|---|---|
| Campo de força | AMBER14ffSB + TIP3P |
| Solvente | TIP3P, 1 nm padding |
| Força iônica | 150 mM NaCl (neutralizado) |
| Temperatura | 300 K (Andersen thermostat) |
| Pressão | 1 bar (NPT — barostato Monte Carlo) |
| Timestep produção | 2 fs |
| Tempo mínimo (validação) | **5 ns** |
| Tempo publicação | **20–100 ns** |
| Equilíbrio | NVT 50 ps → NPT 500 ps + release gradual restraints |
| Frequência DCD | 5.000 steps (~10 ps/frame) |
| Checkpoint | 50.000 steps (~100 ps) |

### 3.3 Plataforma — RTX 3060 CUDA

```
OpenMM 8.2 (env datascience — conda)
Plataforma: CUDA — NVIDIA GeForce RTX 3060 12 GB
Precisão: single (DeterministicForces=false)
Power-cap recomendado: 120 W (nvidia-smi -pl 120, requer admin)
```

**Tempos (sistema ~185k átomos = TcTR+FAD×2+LIG+55k águas+íons):**
- 120 W cap: ~14–15 h / 20 ns
- 161 W (sem cap): ~11 h / 20 ns — **risco de queda de energia (no-break)**

### 3.4 Análise de trajetória (PBC-correto)

**Usar mdtraj v1.11 (env datascience) — NÃO PyMOL para medir distâncias!**
PyMOL não reimageia PBC → distâncias cruas dão ~112 Å (artefato de wrapping).

```python
import mdtraj as md
t = md.load("md_TcTR_inhibitor.dcd", top="minimized_TcTR_inhibitor.pdb")
d = md.compute_distances(t, pairs, periodic=True)  # imagem mínima
```

Script reutilizável: `fig13_3d_kit/analyze_binding_mode.py` (topo+traj, PBC, dithiol auto-detect).

- **RMSD backbone Cα:** ≤ 2,0 Å = estável; run01 atingiu 1,82 ± 0,04 Å (5 ns)
- **RMSD ligante no bolso:** ≤ 2,0 Å = pose mantida; run01 = 0,5 Å médio
- **LIG→FAD (PBC):** run01 (inibidor) = 9,0 Å; Fase B (subversivo) — aguardando
- **LIG→ditiol Cys52/57 (PBC):** run01 = 12,0 Å; ditiol detectado automaticamente como par SG-SG < 2,5 Å mais próximo do FAD

### 3.5 Argumento para revisores (atualizado)

> RMSD do ligante < 2,0 Å durante a produção valida matematicamente a pose — eleva
> o paper de "docking computacional" para "evidência mecanística". O modo inibidor
> (sítio de entrada, RMSD 0,5 Å) e o modo subversivo (si-face FAD) são hipóteses
> complementares que serão discriminadas pela análise comparativa A×B.

### 3.6 Resume pós-queda de energia (protocolo obrigatório)

O script retoma automaticamente do checkpoint se ele existir. Para uma run nova, usar tag distinto:

```bat
:: Fase A — inibidor (retoma do .chk se existir)
python md_simulation.py --fresh --ns 20 --tag TcTR_inhibitor --mode inhibitor

:: Fase B — subversivo (tag distinto obrigatório)
python md_simulation.py --fresh --ns 20 --tag TcTR_subversive --mode subversive

:: Força recomeço do zero (ignora .chk)
python md_simulation.py --fresh --ns 5 --tag TcTR_test --mode inhibitor --no-resume
```

---

## 4. Perfil ADMET

| Propriedade | 2-SCN-NQ | Limite Lipinski | Status |
|---|---|---|---|
| Peso Molecular | 215.23 g/mol | ≤ 500 | ✓ |
| cLogP | 2.41 | ≤ 5 | ✓ |
| Doadores H-bond (HBD) | 0 | ≤ 5 | ✓ |
| Aceitadores H-bond (HBA) | 3 | ≤ 10 | ✓ |
| TPSA | 67.72 Å² | ≤ 140 Å² | ✓ |
| Ligações rotacionáveis | 1 | ≤ 10 | ✓ |
| QED | 0.61 | — | Bom |
| Caco-2 permeabilidade | ~28.5 nm/s | ≥ 20 = alta | ✓ Alta |
| Absorção GI | Alta | — | ✓ |
| BHE penetrante | Não | — | ✓ (desejável para antiparasitário) |
| GHS Classe (in silico) | 3 (Tóxico) | — | Moderada |

Conclusão ADMET: **candidato oral viável** — zero violações Lipinski/Veber.

---

## 5. Toxicidade e Seletividade

### 5.1 Interpretação do Índice de Seletividade (IS)

```
IS = CC50 (célula hospedeira) / IC50 (parasito)
```

| IS | Interpretação | Para o paper |
|---|---|---|
| > 50 | Alta seletividade | Forte argumento para desenvolvimento |
| > 10 | Seletividade favorável | Publicável como argumento de segurança |
| 3–10 | Moderada | Discutir como limitação |
| < 3 | Baixa seletividade | Indica toxicidade não seletiva |

### 5.2 Qual CC50 usar?

| Linhagem | Por que usar | Prioridade |
|---|---|---|
| VERO (rim macaco) | Padrão da literatura antiparasitária | ★ Alta — usar sempre |
| THP-1 (macrófago humano) | O parasito vive DENTRO do macrófago | ★ Alta — argumento clínico |
| J774.A1 (macrófago murino) | Equivalente murino ao THP-1 | Média |
| HepG2 (hepatócito) | Avalia toxicidade hepática oral | Média |

### 5.3 Qual IC50 de parasito usar?

| Forma | Espécie | Relevância clínica | Nota |
|---|---|---|---|
| Tripomastigotas intracelulares | T. cruzi | ★ Mais alta | Forma infectante real em macrófagos |
| Epimastigotas axênicos | T. cruzi | Triagem | Forma vetorial — superestima atividade |
| Amastigotas intracelulares | Leishmania | ★ Mais alta | Forma clínica da LV/LC |
| Promastigotas axênicos | Leishmania | Triagem | Forma vetorial — superestima atividade |

### 5.4 Classificação GHS

| Classe | DL50 oral (rato) | Rótulo | Para o paper |
|---|---|---|---|
| 1–2 | ≤ 50 mg/kg | Fatal | Janela terapêutica estreita |
| **3** | **50–300 mg/kg** | **Tóxico** | IS > 10 justifica progressão com cautela |
| 4 | 300–2000 mg/kg | Nocivo | Perfil favorável |
| 5–6 | > 2000 mg/kg | Pouco/não nocivo | Excelente — forte argumento |

### 5.5 Complexação β-Ciclodextrina

Naftoquinonas têm baixa solubilidade aquosa. A complexação com β-CD pode aumentar 10–100×.

| Parâmetro | Valor estimado | Validação experimental |
|---|---|---|
| Ka (constante de estabilidade) | ~1250 M⁻¹ | DSC, FTIR |
| ΔG de inclusão | ~-4.2 kJ/mol | Solubilidade de fase |
| Aumento de solubilidade | ~12× (a 10 mM β-CD) | HPLC |
| Estequiometria | 1:1 (host:guest) | RMN ¹H |
| Status excipiente β-CD | FDA GRAS — uso oral aprovado | — |

---

## 6. Textos para o Manuscrito

### 6.1 Materiais e Métodos — Docking (EN)

> **Molecular Docking Studies.** The three-dimensional structure of 2-SCN-NQ was generated from its SMILES notation using Open Babel v3.1 (O'Boyle et al., 2011) with the MMFF94 force field (1000 steps conjugate gradient). Crystal structures were retrieved from the Protein Data Bank: TcTR (1BZL), LiTR/proxy *L. amazonensis* (2JK6), hGR (3DJJ), TcCruzain (3GRS), TcFeSOD (1E7W), hMnSOD (1LUV), LmDHODH (3C61), hDHODH (4DVH). Protein preparation included removal of water molecules and buffer ions, while essential cofactors (FAD, FMN, Fe²⁺, Mn²⁺) were retained. Molecular docking was performed using AutoDock Vina v1.2 (Eberhardt et al., 2021) with exhaustiveness = 32. Grid boxes (25 × 25 × 25 Å for reductases and DHODH; 22 × 22 × 22 Å for proteases and SOD) were centered on co-crystallized ligand coordinates. Selectivity was assessed as the ratio Ki(human) / Ki(parasite), where Ki was calculated from the best-ranked pose: Ki = exp(ΔG / RT).

### 6.2 Resultados — Docking (EN, template — preencher com valores reais)

> **Molecular Docking Reveals Preferential Binding to Parasitic Trypanothione Reductase.** 2-SCN-NQ displayed the highest binding affinity for TcTR (ΔG = **X.X** kcal/mol, Ki = **X.X** μM), followed by LiTR (ΔG = **X.X** kcal/mol), while showing markedly lower affinity for the human ortholog hGR (ΔG = **X.X** kcal/mol). The calculated in silico selectivity index (SI = Ki_hGR / Ki_TcTR = **X.X**) is consistent with the experimentally determined selectivity index from cytotoxicity assays. Molecular dynamics simulation over 100 ns confirmed the stability of the TcTR–2-SCN-NQ complex (ligand RMSD = **X.X** ± **X.X** Å), validating the docking pose. These findings support the hypothesis that the antiparasitic activity of 2-SCN-NQ involves preferential inhibition of the trypanothione-dependent antioxidant system, consistent with the ROS overproduction observed experimentally.

### 6.3 Resultados — Docking (PT, template)

> **Docking Molecular Revela Ligação Preferencial à Tripanotiona Redutase Parasitária.** O composto 2-SCN-NQ apresentou maior afinidade de ligação para TcTR (ΔG = **X,X** kcal/mol, Ki = **X,X** μM), seguida por LiTR (ΔG = **X,X** kcal/mol), enquanto exibiu afinidade significativamente menor para o ortólogo humano hGR (ΔG = **X,X** kcal/mol). O índice de seletividade in silico calculado (IS = Ki_hGR / Ki_TcTR = **X,X**) é consistente com o índice de seletividade determinado experimentalmente nos ensaios de citotoxicidade. A simulação de dinâmica molecular confirmou a estabilidade do complexo TcTR–2-SCN-NQ ao longo de 100 ns (RMSD ligante = **X,X** ± **X,X** Å), validando o pose do docking. Estes resultados suportam a hipótese de que a atividade antiparasitária da 2-SCN-NQ envolve inibição preferencial do sistema antioxidante tiol-dependente, consistente com a superprodução de EROs observada experimentalmente.

### 6.4 Discussão — Mecanismo de ação (EN, template)

> The molecular docking analysis provided a mechanistic rationale for the observed selectivity of 2-SCN-NQ. The preferential binding to trypanothione reductase — absent in the mammalian host and essential for thiol-based redox homeostasis in trypanosomatids — aligns with our biochemical findings of increased reactive oxygen species and [oxidative stress markers]. The thiocyanate group at position 2 may act as a soft electrophile interacting with catalytic cysteine residues in TcTR, analogous to other naphthoquinone inhibitors of this enzyme. The lower binding affinity for human GR provides a structural basis for the selectivity index observed in vitro, attributed to the wider substrate-binding cleft of TR that accommodates the spermidine bridge of trypanothione.

---

## 7. Narrativa Mecanística Central (atualizada 2026-05-29)

### 7.1 Dois mecanismos possíveis — não mutuamente exclusivos

**A. INIBIDOR DE TR (hipótese A, suportada por docking+MD atual)**

O 2-SCN-NQ ocupa **estavelmente** (RMSD 0,5 Å, 5 ns) a **fenda hidrofóbica de entrada** do sítio ativo da TcTR — o mesmo bolso onde se ligam inibidores clássicos como mepacrina e fenotiazinas (sítio tipo "Z"). Distância ao FAD ~9 Å, ao ditiol ~12 Å.

- Mecanismo: bloqueia o acesso da tripanotiona ao sítio catalítico → inibição competitiva
- Seletividade vs hGR: estrutural (canal TR mais largo; Ala34 onde hGR tem Arg)
- Análogos documentados: mepacrina, clorpromazina, tioridazina, série NQ-2-substituída

**B. SUBSTRATO SUBVERSIVO (hipótese B, em teste — Fase B MD)**

O núcleo quinônico empilha na **si-face da isoaloxazina do FAD** (pose construída a 3,60 Å). A quinona aceita elétrons do FADH₂ → semiquinona → ciclo redox com O₂ → superóxido → estresse oxidativo.

- Mecanismo: "sequestro" do fluxo eletrônico da TR → inversão: defesa → ataque
- Dados experimentais que suportam: DHE↑ (55% DHE+ a 2×IC₅₀), colapso ΔΨm, −80% respiração
- Status: hipótese B aguarda resultado da MD da si-face do FAD (tag `TcTR_subversive`)

**C. DUAL (drug design em andamento)**

A molécula atual pode ser pequena demais para tocar os dois sítios simultaneamente.
O pipeline `dual_design_pipeline.py` propõe 5 análogos maiores com cauda de poliamina/
acridina/fenotiazina que alcançam ambos os sub-sítios.

### 7.2 Dados experimentais reais (2-SCN-NQ, projeto 2015–2021)

| Parâmetro | Valor | Método | Relevância |
|---|---|---|---|
| IC₅₀ T. cruzi epimastigotas | 5,5 µM | MTT 24h, Dm28c | Triagem |
| IC₅₀ T. cruzi tripomastigotas | **0,3 µM** | Mobilidade 24h, cepa Y | Forma infectante ★ |
| IC₅₀ T. cruzi amastigotas | 4,72 µM | Vero infectadas 48h | Forma clínica ★ |
| IC₅₀ L. amazonensis promastigotas | 5,64 µM (24h) | MTT, L. amazonensis | Triagem |
| CC₅₀ Vero | >50 µM | MTT 48h | Menos sensível |
| CC₅₀ RAW264.7 (macrófago murino) | 25 µM | MTT 48h | Referência toxicidade |
| CC₅₀ HeLa | 10,36 µM | MTT | Mais sensível |
| IS tripomastigotas/Vero | **333** | — | Excelente |
| IS epimastigotas/Vero | 18 | — | Favorável |
| IS Leishmania/Vero | 17,7 | — | Favorável |
| Benznidazol IC₅₀ (controle) | 16,3 µM | MTT 24h, cepa Y | Mesmo protocolo |

### 7.3 Argumento central (paper ACS Inf Dis)

> O 2-SCN-NQ liga-se **estavelmente** à TcTR no sítio de entrada (inibidor) **e/ou** na
> si-face do FAD (substrato subversivo), gerando um duplo bloqueio: impede o acesso da
> tripanotiona **e** gera superóxido que o parasito não consegue neutralizar (TR bloqueada).
> A seletividade estrutural (Ala34/TR vs Arg/hGR) explica a janela terapêutica experimental
> (IS 18–333). Dados: ROS experimental (DHE↑) + docking (TcTR > hGR) + MD estável + IS ✓ ADMET ✓.

---

## 8. Checklist para Submissão

### Pipeline computacional
- [ ] exhaustiveness ≥ 32 (para reprodutibilidade)
- [ ] Grid box centrado corretamente (conferir co-cristal pelo QC do pipeline)
- [ ] Cofatores FAD/FMN/metais retidos nas proteínas
- [ ] Comparação parasito vs humano em ≥ 1 par de ortólogos
- [ ] MD confirmou estabilidade (RMSD ≤ 2.0 Å em ≥ 80% da trajetória de produção)

### Dados experimentais (aba Toxicidade — auto-carregados p/ 2-SCN-NQ)
- [x] IC50 tripomastigotas T. cruzi = **0,3 µM** (mobilidade, cepa Y)
- [x] IC50 epimastigotas T. cruzi = **5,5 µM** (MTT 24h, Dm28c)
- [x] IC50 amastigotas intracelulares T. cruzi = **4,72 µM** (Vero 48h)
- [x] IC50 L. amazonensis promastigotas = **5,64 µM** (24h)
- [x] CC50 VERO = **>50 µM** (sem citotoxicidade significativa <10 µM)
- [x] CC50 RAW264.7 (J774) = **25 µM** (macrófago murino)
- [x] CC50 HeLa = **10,36 µM**
- [x] Benznidazol IC50 = **16,3 µM** (controle, mesmo protocolo)
- [x] IS tripomastigotas/Vero = **333** | IS epimastigotas/Vero = **18** | IS Leish/Vero = **17,7**
- [ ] DL50 camundongo oral (mg/kg) — não determinada formalmente

### Figuras para o manuscrito
- [ ] Fig. 1: Barras comparativas ΔG — todos os alvos (PNG 300 dpi + TIFF para ACS)
- [ ] Fig. 2: Índices de seletividade (SI) — parasito vs humano
- [ ] Fig. 3: Heatmap ΔG + Ki — todos os alvos
- [ ] Fig. 4: Painel RMSD/RMSF da MD (100 ns, alvo principal)
- [ ] Fig. 5: Perfil ADMET (BOILED-Egg ou dashboard figure)
- [ ] Suplementar: Poses PyMOL dos 3 melhores alvos

### Texto do manuscrito
- [ ] Métodos: citar versões de todos os softwares
- [ ] Métodos: descrever critério de manutenção de cofatores
- [ ] Métodos: reportar RMSD de re-docking (validação do protocolo) < 2.0 Å
- [ ] Resultados: conectar IS in silico com IS experimental
- [ ] Discussão: mecanismo dual ROS + inibição TR
- [ ] Dados suplementares: coordenadas dos grid boxes, todas as poses top-ranked

### Formatos para ACS
- Figuras: TIFF 300 dpi (gerado automaticamente pelo pipeline)
- Tabela: CSV → importar no Word/Excel
- SI: RMSD dos re-dockings, grid box coordinates, compound SMILES

---

## 9. Referências Científicas dos Alvos

| Alvo | PDB | Referência chave |
|---|---|---|
| TcTR (1BZL) | Melby et al. (2003) *J. Biol. Chem.* 278:39694 |
| LiTR proxy (2JK6) | Baiocco et al. (2009) *Structure* 17:1547 |
| TcCruzain (3GRS) | Sajid & McKerrow (2002) *Mol. Biochem. Parasitol.* 120:1 |
| TcFeSOD (1E7W) | Wilkinson et al. (2000) *J. Biol. Chem.* 275:15 |
| LmDHODH (3C61) | Rodrigues et al. (2021) *Molecules* 26:6002 |
| hGR (3DJJ) | Berkholz et al. (2008) *J. Mol. Biol.* 376:1444 |
| hCatL (2XU3) | Nagler et al. (1997) *Biochemistry* 36:12608 |
| hMnSOD (1LUV) | Borgstahl et al. (1992) *Biochemistry* 31:4317 |
| hDHODH (4DVH) | Baumgartner et al. (1996) *J. Mol. Biol.* 259:14 |
| Proxy *L. amazonensis* | Camargo et al. (2023) *Sci. Rep.* 13:6887 |
