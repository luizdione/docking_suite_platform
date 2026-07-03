# Deploy do NQ-Docking Viewer na Hostinger

Versão **hospedada** e **somente leitura** da plataforma docking. Publica **exatamente os
valores, tabelas e figuras do manuscrito ACS submetido** do 2-SCN-NQ (aba **Manuscrito** com
o PDF preliminar embutido) + uma **calculadora ADMET que roda no navegador**. É um site 100%
**estático** — funciona em qualquer plano Hostinger, sem Python/Node/GPU no servidor.

> **Onde está publicado hoje:** `https://ifrj-crj-geneticamolecular.online/docking/`
> (subpasta do domínio principal — os caminhos do viewer são relativos e funcionam nesse
> layout). Se preferir um subdomínio `docking.ifrj-crj-geneticamolecular.online`, crie-o no
> hPanel → Subdomínios e envie o mesmo conteúdo para a pasta raiz dele.

## 0. Gerar o pacote (na sua máquina)

```powershell
# conda datascience; a partir da raiz do repo docking/
python export_web.py --out dist_web
# ativar a senha (Basic Auth) — a senha vira hash; nunca é gravada em texto:
python enable_basic_auth.py --out dist_web --user darwin --password "SUA_SENHA"
```

O exportador **lê os dados científicos do manuscrito** (`2-SNC-NQ/PAPER/…`, via `NQ_PAPER_DIR`)
e copia as **figuras oficiais** (fig10–13, estrutura) e o **PDF de submissão**. O `dist_web/`
resultante contém: `index.html`, `lib/`, `.htaccess`, `.htpasswd` (se ativou a senha),
`config.json`, `data/`, `figures/`, `assets/` (logo, pipeline, PDF) e `docs/`.

## 1. Enviar os arquivos

Envie **o conteúdo** de `dist_web/` (não uma pasta por fora) para o document root — hoje
`.../public_html/docking/`. Pelo **Gerenciador de Arquivos** do hPanel: compacte o conteúdo num
`.zip`, faça upload e **extraia lá dentro**; garanta que `index.html` fique na raiz da pasta.
Inclua os arquivos ocultos `.htaccess` e `.htpasswd`.

## 2. Senha (Basic Auth) — usuário `darwin`

Duas formas:

- **(A) Arquivos gerados** (feito por `enable_basic_auth.py`): o `.htaccess` já vem com
  `AuthType Basic` + `AuthUserFile` apontando para o caminho absoluto do `.htpasswd`.
  Se der **erro 401/500**, o caminho do `AuthUserFile` está errado: acesse
  `.../docking/_auth_path.php` **uma vez**, copie o caminho que ele imprime para o
  `AuthUserFile` do `.htaccess`, salve e **apague o `_auth_path.php`**.
- **(B) GUI do hPanel** (infalível, sem mexer em caminho): hPanel → **Avançado** →
  **Proteção por senha de diretórios** → selecione a pasta `docking` → crie o usuário `darwin`
  com a senha. Se usar esta via, pode remover o `.htpasswd`/bloco Basic Auth do `.htaccess`.

## 3. Atualizar

Rode `python export_web.py --out dist_web` (+ `enable_basic_auth.py`) de novo e reenvie
`data/`, `figures/`, `assets/` e `config.json`. O `index.html`/`lib/` raramente mudam.

## Fidelidade ao paper

Os dados de **docking** (ΔG/Ki/PDB por alvo), **MD** (distâncias ligante–FAD n=3, contatos,
estatística) e **ADMET** vêm do manuscrito submetido (`2scnnq_manuscript_EN.tex`), não dos dados
de trabalho. As figuras são as **oficiais** do paper (`PAPER/figures/`). A aba **Manuscrito**
mostra o PDF de submissão embutido, para conferência direta.

## O que roda ONLINE (sem servidor)

| Aba | No Hostinger |
|---|---|
| Visão geral · Docking · Dinâmica Molecular · Manuscrito · Subprojetos · Ajuda | leitura (dados/figuras do paper) |
| **ADMET · Compostos** | **leitura + calculadora client-side (online)** |

Camada B opcional (SMILES → descritores): coloque `lib/RDKit_minimal.js` e `.wasm` para habilitar.
