# =============================================================================
# docking_suite_platform — imagem para piloto de monetizacao (RunPod)
# Ver .project-memory/DECISIONS.md [D-0010] (OneDrive\Projetos\docking) para o
# racional de negocio. Esta imagem roda O SERVICO, nao e distribuida ao cliente.
#
# Base CUDA (nao conda-forge): o requirements.txt do projeto e 100% pip-instalavel
# (openmm>=8.4 publica wheel PyPI com suporte CUDA desde a serie 8.x; RDKit e
# openbabel-wheel idem). Isso replica o metodo "ALTERNATIVO (pip)" ja documentado
# no proprio requirements.txt — nao precisa de conda aqui.
# =============================================================================
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONUTF8=1 \
    PIP_NO_CACHE_DIR=1 \
    PATH="/opt/venv/bin:$PATH"

# --- Python 3.11 + libs de sistema exigidas por RDKit/OpenBabel/OpenMM ---------
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3-pip \
        wget unzip ca-certificates \
        libxrender1 libxext6 libsm6 build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

WORKDIR /app

# --- Criação do Ambiente Virtual Isolado para mitigar erros de dependência -----
RUN python3 -m venv /opt/venv

# --- Dependencias Python (camada cacheavel) ------------------------------------
COPY requirements.txt .

# Utiliza o resolver antigo emparelhado com o ambiente isolado do venv
RUN pip install --upgrade pip setuptools wheel && \
    pip install --use-deprecated=legacy-resolver -r requirements.txt

# --- AutoDock Vina (nao instalavel via pip — baixado do release oficial) ------
# Binario Linux oficial do AutoDock Vina (ccsb-scripps). O valor anterior baixava
# a HOME do github.com como "vina" -- a imagem nunca rodaria docking. --version
# no build FALHA a imagem se o binario nao for executavel (sanidade fail-loud).
ARG VINA_VERSION=1.2.5
RUN wget -qL "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
        -O /usr/local/bin/vina && \
    chmod +x /usr/local/bin/vina && \
    /usr/local/bin/vina --version

# --- Codigo da plataforma -------------------------------------------------------
COPY . .

# Dados: o build LOCAL embute docking_project/ via COPY . . (fica em disco mesmo
# fora do git). O codigo, porem, procura em E:\data\docking por default (path
# Windows que NAO existe aqui) -- estas envs redirecionam para os dados embutidos.
# Para rodar com dados montados (ex.: volume no RunPod), sobreponha em runtime:
#   docker run -e DOCKING_DATA_DIR=/data -v /host/dados:/data ...
ENV DOCKING_DATA_DIR=/app/docking_project \
    NQ_DATA_DIR=/app/docking_project

# Porta do dashboard Flask (dashboard_v2.py aceita --host/--port)
EXPOSE 5001

# ENTRYPOINT dual:
#   - modo servico (piloto): 1 job de cada vez via runner.py (ja tem CLI de
#     proveniencia: plan/run/status/resume) — e o que o RunPod deve chamar por job.
#   - modo dashboard: para inspecao manual/relatorio, nao para expor ao cliente.
# Ex.: docker run --gpus all IMG python3 runner.py run --campaign SUB01 --kind docking
ENTRYPOINT ["python3"]
CMD ["dashboard_v2.py", "--host", "0.0.0.0", "--port", "5001"]
