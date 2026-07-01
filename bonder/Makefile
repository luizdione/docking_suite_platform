# bonder - atalhos de desenvolvimento e operacao
.PHONY: help install test up up-gpu down dbsync baits qc lint

help:
	@echo "bonder - alvos disponiveis:"
	@echo "  make install    instala o pacote + deps de dev (pip -e .[dev,web,ms])"
	@echo "  make test       roda a suite de testes (nucleo CPU, sem GPU/R)"
	@echo "  make up         sobe a plataforma (CPU) em http://localhost:8080"
	@echo "  make up-gpu     sobe a plataforma com overlay de GPU (RTX 3060)"
	@echo "  make down       derruba os conteineres"
	@echo "  make dbsync     baixa/verifica os bancos murinos (bonder db sync)"
	@echo "  make baits      valida iscas e monta o FASTA de busca"

install:
	pip install -e ".[dev,web,ms,annotate]"

test:
	python -m pytest -q

up:
	docker compose -f docker/docker-compose.yml up --build

up-gpu:
	docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build

down:
	docker compose -f docker/docker-compose.yml down

dbsync:
	python -m bonder.cli db sync

baits:
	python -m bonder.cli baits check
