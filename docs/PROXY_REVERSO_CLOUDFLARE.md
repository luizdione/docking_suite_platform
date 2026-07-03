# Proxy Reverso — expor o dashboard completo via Cloudflare Tunnel

Este guia descreve como publicar o **dashboard completo** (`dashboard_v2.py`, `localhost:5001`,
com docking/MD/GPU) na internet, de forma segura, a partir da sua máquina com **IP dinâmico**
(fibra VIVO, sem IP fixo e sem abrir portas no roteador), usando **Cloudflare Tunnel**
(`cloudflared`) e o domínio `ifrj-crj-geneticamolecular.online`.

> **Quando usar cada versão**
>
> | Necessidade | Solução |
> |---|---|
> | Publicar **resultados** (leitura), sempre no ar, sem depender da sua máquina | **Hostinger estático** (`export_web.py` → `docking.ifrj-crj-geneticamolecular.online`) — ver `web/viewer/README_DEPLOY.md` |
> | Acessar remotamente o **app completo** (rodar docking/MD/GPU, gerar figuras) | **Cloudflare Tunnel** (este guia) — só funciona com a máquina ligada |
>
> As duas abordagens coexistem: `docking.` (estático, leitura) e, por exemplo, `app-docking.`
> (túnel, completo). O túnel **não** contorna a GPU — respeite o power-cap e o gate de GPU.

## Por que Cloudflare Tunnel (e não port forwarding)

- **IP dinâmico não é problema**: o `cloudflared` abre uma conexão **de saída** para a
  Cloudflare; não há porta de entrada aberta no roteador, nem necessidade de DDNS/IP fixo.
- **HTTPS automático** (certificado gerenciado pela Cloudflare).
- **Autenticação na borda** (Cloudflare Access): exige login antes de chegar ao Flask —
  essencial, pois o dashboard **não** tem autenticação própria e dispara jobs pesados.

## Pré-requisitos

1. O domínio `ifrj-crj-geneticamolecular.online` com **nameservers apontados para a Cloudflare**
   (plano gratuito serve). No hPanel da Hostinger, troque os NS pelos que a Cloudflare indicar
   ao adicionar o site (Painel Cloudflare → *Add a site*). A propagação leva de minutos a horas.
2. `cloudflared` instalado na máquina (Windows):
   ```powershell
   winget install --id Cloudflare.cloudflared
   # ou baixe cloudflared.exe de https://github.com/cloudflare/cloudflared/releases
   ```

## Passo a passo

### 1. Autenticar o cloudflared na sua conta Cloudflare
```powershell
cloudflared tunnel login
```
Abre o navegador; escolha a zona `ifrj-crj-geneticamolecular.online` e autorize.

### 2. Criar o túnel e a rota DNS
```powershell
cloudflared tunnel create docking-app
# cria um Tunnel ID e um arquivo de credenciais em %USERPROFILE%\.cloudflared\<UUID>.json
cloudflared tunnel route dns docking-app app-docking.ifrj-crj-geneticamolecular.online
```

### 3. Arquivo de configuração
Crie `%USERPROFILE%\.cloudflared\config.yml`:
```yaml
tunnel: docking-app
credentials-file: C:\Users\luizd\.cloudflared\<UUID>.json
ingress:
  - hostname: app-docking.ifrj-crj-geneticamolecular.online
    service: http://localhost:5001
  - service: http_status:404
```

### 4. Subir o túnel (com o dashboard já rodando em :5001)
```powershell
python dashboard_v2.py           # terminal 1 (idealmente como Administrador, p/ power-cap)
cloudflared tunnel run docking-app   # terminal 2
```
Acesse: `https://app-docking.ifrj-crj-geneticamolecular.online`.

### 5. (Fortemente recomendado) Proteger com Cloudflare Access
Painel Cloudflare → **Zero Trust** → **Access** → **Applications** → *Add* → *Self-hosted*:
- Domínio: `app-docking.ifrj-crj-geneticamolecular.online`.
- Política: *Allow* apenas o(s) e-mail(s) autorizado(s) (ex.: `luiz.melo@ifrj.edu.br`), via
  código enviado por e-mail (One-time PIN) ou provedor de identidade.

Sem esse passo, **qualquer pessoa com a URL** teria acesso ao dashboard e poderia disparar MD/GPU.

## Rodar como serviço (opcional, para ficar sempre no ar)
```powershell
cloudflared service install       # PowerShell como Administrador
```
Instala o túnel como serviço do Windows (sobe no boot). Lembre que a máquina precisa estar ligada
e o dashboard rodando para o app responder.

## Avisos de segurança e GPU

- **Autenticação é obrigatória** (Cloudflare Access): o Flask não tem login e as rotas `/run_*`
  disparam docking/MD. Nunca exponha o túnel sem Access.
- **GPU compartilhada / no-break**: acesso remoto não muda as regras locais — não dispare MD com o
  LILO/Ollama carregado; aplique o power-cap (`nvidia-smi -pl 120`) e use o `run_md_safe.ps1`.
- **Segredos**: `ANTHROPIC_API_KEY` e afins continuam apenas em variáveis de ambiente locais;
  nunca no repositório nem no bundle publicado.
- Prefira publicar **resultados** pela via estática (Hostinger) e reservar o túnel para quando
  precisar operar o app remotamente.
