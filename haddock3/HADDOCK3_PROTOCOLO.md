# HADDOCK3 — docking guiado por informação do 2-SCN-NQ na TcTR

Docking *information-driven* como método **ortogonal/consenso** ao AutoDock Vina,
ancorado no sítio do substrato (trypanothione) da *T. cruzi* trypanothione
reductase (estrutura 1BZL). Responde à crítica da banca de incorporar HADDOCK.

## Por que HADDOCK aqui (resumo de defesa)

O Vina é apropriado para molécula pequena em bolso, mas é *ab initio* (cego). O
HADDOCK usa **restrições de interação ambíguas (AIRs)** a partir do conhecimento
do sítio ativo — aqui, os resíduos que coordenam a glutationa/trypanothione
co-cristalizada. Reportar **consenso Vina + HADDOCK** (pose comum + score
HADDOCK + RMSD entre os dois) é mais robusto do que um método só.

## Arquivos gerados (nesta pasta)

| Arquivo | Conteúdo |
|---|---|
| `TcTR_receptor.pdb` | Receptor (dímero 1BZL, cadeias A/B, 7459 átomos, segid = cadeia) |
| `ligand_2SCN-NQ.pdb` | Ligante (pose 1 do Vina; chain C, resíduo `LIG`; nomes de átomo únicos) |
| `ambig.tbl` | AIRs: CYS53, CYS58, HIS461 (+ bolso) ↔ ligante |
| `active_site_residues.json` | Resíduos do sítio (< 4,5 Å da glutationa co-cristalizada) |
| `run-TcTR-2SCN-NQ.cfg` | Workflow HADDOCK3 (topoaa → rigidbody → flexref → emref → clustfcc) |

> Nota de numeração: no 1BZL as cisteínas catalíticas são **CYS53/CYS58**
> (= Cys52/Cys57 da numeração canônica; o `docking_v2.py` usa 52/57, +1 de
> offset). HIS461 é a histidina catalítica do parceiro do dímero. As AIRs já
> usam a numeração correta do 1BZL.

## 1. Instalar o HADDOCK3 (Linux/WSL)

HADDOCK3 é Linux/Mac. No Windows, use **WSL2** (Ubuntu). Requer CNS (licença
acadêmica gratuita).

```bash
# em WSL/Ubuntu
python3 -m venv haddock3-env && source haddock3-env/bin/activate
pip install haddock3
# CNS: solicitar binario academico em https://cns-online.org/ e exportar:
export HADDOCK3_CNS_EXEC=/caminho/para/cns_solve
haddock3 --version
```

## 2. ⚠️ Passo crítico — topologia/parâmetros CNS do ligante

O `LIG` é um resíduo não-padrão; o HADDOCK precisa de `lig.top` (topologia) e
`lig.param` (parâmetros) em formato CNS. Use o `2-SCN-NQ.mol2` que vocês já têm
(`docking_project/ligand/2-SCN-NQ/2-SCN-NQ.mol2`). Rotas:

- **Recomendada — tutorial oficial de ligantes do HADDOCK** (gera .top/.param a
  partir do mol2/PDB):
  https://www.bonvinlab.org/education/HADDOCK3/HADDOCK3-protein-ligand/
- **PRODRG / ATB** (Automated Topology Builder, https://atb.uq.edu.au): submeta
  o mol2, baixe a topologia e converta para CNS.
- **acpype** (GAFF2/AM1-BCC) → converter para CNS com `ccpn`/scripts da comunidade.

Garanta que os **nomes de átomo** no `lig.top` batem com os do
`ligand_2SCN-NQ.pdb` (C1, C2, …, O1, O2, S1, N1). Coloque `lig.top`/`lig.param`
nesta pasta.

## 3. Rodar

```bash
cd haddock3
haddock3 run-TcTR-2SCN-NQ.cfg
```

## 4. Analisar / reportar (consenso com o Vina)

- Saída em `run-TcTR-2SCN-NQ/` — o módulo `caprieval` final dá o **HADDOCK
  score** por cluster e RMSD intra-cluster; o `analysis/` traz plots.
- **RMSD Vina↔HADDOCK:** alinhe o receptor e calcule o RMSD heavy-atom entre a
  melhor pose HADDOCK (cluster 1) e a pose 1 do Vina
  (`docking_project/docking/result_TcTR.pdbqt`). RMSD baixo = **consenso** entre
  os dois métodos → resultado defensável. Pode reaproveitar
  `docking_pose_rmsd.py` (mesma métrica) para esse cálculo.

## 5. O que reportar no paper/banca

1. Vina (ab initio): afinidade + clustering de poses por RMSD (`docking_pose_rmsd_TcTR.png`).
2. HADDOCK (guiado por AIRs do sítio catalítico): HADDOCK score + cluster top.
3. **Consenso:** RMSD entre as melhores poses dos dois métodos; contatos comuns
   com CYS53/CYS58/HIS461.
4. Validação do protocolo: redocking RMSD < 2 Å (`docking_rmsd_validation.py`).
