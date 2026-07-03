# Deploy do NQ-Docking Viewer na Hostinger (docking.ifrj-crj-geneticamolecular.online)

Versão **hospedada** e **leve** da plataforma docking: **somente leitura** dos
resultados já gerados (docking, binding mode, MD/RMSD, ADMET, compostos,
subprojetos) + **uma calculadora de drug-likeness que roda no próprio navegador**.
A análise pesada (docking com Vina, dinâmica molecular OpenMM/GPU, geração de
figuras) continua rodando **local** (dashboard em `localhost:5001`); aqui só
publicamos os resultados. É um site 100% **estático** — funciona em qualquer plano
Hostinger (inclusive o compartilhado), sem Python/Node/GPU no servidor.

## 0. Gerar o pacote (na sua máquina)

```powershell
# conda datascience; a partir da raiz do repo docking/
python export_web.py --out dist_web
# opcional: apenas alguns alvos de MD ou outro composto
python export_web.py --out dist_web --compound 2-SCN-NQ
```

O `dist_web/` gerado já contém: `index.html`, `lib/`, `.htaccess`, `config.json`,
`data/` (compostos, ADMET, binding_mode, md/*, subprojetos), `figures/`,
`assets/pipeline_SUB01.svg` e `docs/` (com `index.json` + os .md da Ajuda).

## 1. Criar o subdomínio (hPanel)

1. hPanel → **Domínios** → **Subdomínios**.
2. Crie `docking` para `ifrj-crj-geneticamolecular.online`.
3. Anote a **pasta raiz (document root)**, ex.:
   `domains/ifrj-crj-geneticamolecular.online/public_html/docking`
   (a Hostinger cria o registro DNS automaticamente; pode levar alguns minutos).

## 2. Enviar os arquivos

- **Gerenciador de Arquivos (hPanel)**: compacte o **conteúdo** de `dist_web/` num
  `.zip`, envie para a pasta raiz do subdomínio e **extraia lá dentro**.
  Garanta que `index.html` fique na RAIZ do subdomínio (não dentro de uma subpasta).
- ou **FTP/SFTP** (FileZilla/WinSCP): copie tudo de `dist_web/` para a pasta raiz.

Confira: `https://docking.ifrj-crj-geneticamolecular.online/` deve abrir o viewer.

## 3. Proteger por senha (a "área interna")

**(A) Recomendado — GUI:** hPanel → **Avançado** →
**Proteção por senha de diretórios** → selecione a pasta do subdomínio →
crie usuário e senha. Pronto (a Hostinger cuida do `.htaccess`/`.htpasswd`).

**(B) Manual:** edite o `.htaccess` do bundle, descomente o bloco *Basic Auth* e
aponte `AuthUserFile` para o caminho ABSOLUTO do seu `.htpasswd` no servidor.

> O `.htaccess` incluído já vem **seguro por padrão** (não ativa a senha sozinho,
> para evitar erro 500); ele só define MIME correto (inclui `.wasm`), desliga a
> listagem de diretórios e cacheia os assets. A senha é ativada por (A) ou (B).

## 4. Atualizar os resultados (nas próximas corridas)

Rode `python export_web.py --out dist_web` de novo e reenvie **apenas** as pastas
`data/`, `figures/` e o `config.json` (o `index.html`/`lib/` raramente mudam).

## O que roda ONLINE (sem servidor) vs. só LOCAL

| Módulo / aba              | No Hostinger (estático)                        | Observação |
|---------------------------|------------------------------------------------|------------|
| Visão geral / Seletividade| **leitura** (resultados exportados)            | — |
| Docking & Binding         | **leitura** (binding_mode.json)                | docking com Vina roda local |
| MD (RMSD/RMSF)            | **leitura** (resumos + figuras)                | MD OpenMM exige **GPU** → só local |
| **ADMET**                | **leitura + calculadora client-side (online)** | regras QSAR reimplementadas em JS |
| **Compostos**            | **leitura + descritores → calculadora**        | validação de SMILES via RDKit.js opcional |
| Subprojetos / Ajuda       | **leitura** (SVG + docs .md)                    | — |

### Camada B opcional — SMILES → descritores no navegador (RDKit.js/WASM)

A calculadora ADMET já roda online a partir dos descritores (MW/logP/HBD/HBA/TPSA/RB).
Para aceitar **um SMILES novo** e derivar os descritores automaticamente (+ desenho 2D),
coloque os dois arquivos do [RDKit MinimalLib](https://www.npmjs.com/package/@rdkit/rdkit)
em `lib/`:

```
lib/RDKit_minimal.js
lib/RDKit_minimal.wasm
```

O viewer detecta e habilita o botão **"SMILES → descritores"** automaticamente.
Sem esses arquivos, a calculadora continua funcionando com entrada manual de descritores.
