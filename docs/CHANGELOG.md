# Changelog — NQ-Docking Suite

Formato: [versão] — data. Projeto: docking + MD da TcTR (T. cruzi) com 2-SCN-NQ.

## [2.1.0] — 2026-05-22

Resposta às críticas de banca (RMSD + HADDOCK) e correção definitiva do force-field do FAD.

### Adicionado
- **Validação de docking por RMSD** (`docking_pose_rmsd.py`): clustering das poses
  do Vina por RMSD pareado + figura por alvo (`docking_project/figures/docking_pose_rmsd_*.png`).
- **Validação do protocolo por redocking** (`docking_rmsd_validation.py`): re-docka a
  glutationa co-cristalizada (GCG) do 1BZL e mede RMSD vs. cristal (torna real o
  "RMSD < 2 Å" do texto de Métodos).
- **HADDOCK3** (`haddock3/`): docking guiado por informação como método ortogonal —
  receptor (dímero 1BZL), ligante 2-SCN-NQ, AIRs ancoradas em CYS53/CYS58/HIS461,
  workflow `run-TcTR-2SCN-NQ.cfg` e protocolo `HADDOCK3_PROTOCOLO.md`.
- **Plataforma (dashboard) v2.1**: nova aba "RMSD & HADDOCK" + rota `/fig_rmsd/<fig>`.
- Versionamento: `.gitignore`, `VERSION`, `CHANGELOG.md`.

### Corrigido (force-field do FAD — destrava a MD parada desde 14/05)
- `regen_fad_xml.py`: regenera o force-field do FAD com elementos/massas/LJ/cargas/
  bonds/angles/**torsions** corretos (rota SMIRNOFF + ordens de ligação via PDB CCD,
  determinístico — sem perceção de valência que falhava no fósforo).
- `cofactor_FAD_TcTR_v3.pdb`: coluna `element` corrigida (33 átomos).
- `md_clean.py`: passa a preferir o par v3 (XML+PDB) quando presente, com fallback
  consistente para o antigo.
- `auto_fix/fix_fad_xml.py`: corrigido bug que injetava 91 bonds duplicadas
  (regex `type1/type2`); marcado como stopgap.

### Diagnóstico (documentado)
- `BRIEFING_2026-05-21_md_module.md`, `FIX_APPLIED_2026-05-21.md`: o XML antigo do FAD
  tinha 33/86 átomos com elemento/LJ errados (scramble) e `PeriodicTorsionForce` vazio.

## [2.0.0] — anterior
- Pipeline docking (Vina) + ADMET + MD (OpenMM) + ciclodextrina + toxicidade + dashboard Flask.
