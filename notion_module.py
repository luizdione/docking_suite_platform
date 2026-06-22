#!/usr/bin/env python3
"""
notion_module.py — Integração NQ-Docking Suite ↔ Notion
========================================================
Funções para:
  • setup_notion_databases()   — cria databases de Experimentos, Compostos e Referências
  • log_experiment()           — registra análise (ADMET, docking, tox) no Notion
  • export_report()            — exporta relatório completo como página Notion
  • sync_compound()            — sincroniza/atualiza composto no database
  • import_compounds()         — importa compostos do Notion para o projeto
  • import_references()        — importa referências bibliográficas do Notion

Configuração:
  Requer NOTION_TOKEN no ambiente ou em docking_project/notion_config.json
  A página-pai padrão é a página "docking" (ID: 34a34d23-b07f-8140-9002-d06ddda48305)

Uso:
  import notion_module
  notion_module.setup_notion_databases()
  notion_module.log_experiment(session_tag, compound_name, "admet", admet_data)
"""

import os
import sys
import json
import datetime
from pathlib import Path
from typing import Optional

# ── Configuração ────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
CONFIG_FILE  = SCRIPT_DIR / "docking_project" / "notion_config.json"
NOTION_PAGE_DOCKING = "34a34d23-b07f-8140-9002-d06ddda48305"

# IDs dos databases criados (populados por setup_notion_databases)
DB_EXPERIMENTS  = None
DB_COMPOUNDS    = None
DB_REFERENCES   = None


# ── Helpers de configuração ─────────────────────────────────────────────────

def _get_token() -> Optional[str]:
    """Retorna o token Notion da variável de ambiente ou do config JSON."""
    token = os.environ.get("NOTION_TOKEN")
    if token:
        return token
    if CONFIG_FILE.exists():
        try:
            cfg = json.loads(CONFIG_FILE.read_text())
            return cfg.get("token")
        except Exception:
            pass
    return None


def _get_db_ids() -> dict:
    """Carrega IDs dos databases do arquivo de configuração."""
    if CONFIG_FILE.exists():
        try:
            cfg = json.loads(CONFIG_FILE.read_text())
            return cfg.get("databases", {})
        except Exception:
            pass
    return {}


def _save_db_ids(db_ids: dict):
    """Salva IDs dos databases no arquivo de configuração."""
    CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    cfg = {}
    if CONFIG_FILE.exists():
        try:
            cfg = json.loads(CONFIG_FILE.read_text())
        except Exception:
            pass
    cfg.setdefault("databases", {}).update(db_ids)
    CONFIG_FILE.write_text(json.dumps(cfg, indent=2, ensure_ascii=False))
    print(f"[Notion] Config salvo em {CONFIG_FILE}")


def _notion_headers(token: str) -> dict:
    return {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
        "Notion-Version": "2022-06-28",
    }


def _notion_request(method: str, endpoint: str, payload: dict = None) -> dict:
    """Faz requisição HTTP à API Notion."""
    import urllib.request
    import urllib.error

    token = _get_token()
    if not token:
        raise ValueError(
            "Token Notion não encontrado. Defina NOTION_TOKEN no ambiente "
            "ou adicione {'token': 'secret_...'} em docking_project/notion_config.json"
        )

    url = f"https://api.notion.com/v1/{endpoint}"
    headers = _notion_headers(token)
    data = json.dumps(payload).encode("utf-8") if payload else None

    req = urllib.request.Request(url, data=data, headers=headers, method=method.upper())
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8")
        raise RuntimeError(f"Notion API {method} /{endpoint} → {e.code}: {body}")


def _text_block(content: str, block_type: str = "paragraph") -> dict:
    """Cria um bloco de texto Notion."""
    return {
        "object": "block",
        "type": block_type,
        block_type: {
            "rich_text": [{"type": "text", "text": {"content": content[:2000]}}]
        }
    }


def _heading_block(content: str, level: int = 2) -> dict:
    t = f"heading_{level}"
    return {
        "object": "block",
        "type": t,
        t: {"rich_text": [{"type": "text", "text": {"content": content}}]}
    }


def _divider_block() -> dict:
    return {"object": "block", "type": "divider", "divider": {}}


def _callout_block(content: str, emoji: str = "🧪") -> dict:
    return {
        "object": "block",
        "type": "callout",
        "callout": {
            "rich_text": [{"type": "text", "text": {"content": content[:2000]}}],
            "icon": {"type": "emoji", "emoji": emoji}
        }
    }


def _code_block(content: str, language: str = "plain text") -> dict:
    return {
        "object": "block",
        "type": "code",
        "code": {
            "rich_text": [{"type": "text", "text": {"content": content[:2000]}}],
            "language": language
        }
    }


def _table_row(cells: list) -> dict:
    return {
        "object": "block",
        "type": "table_row",
        "table_row": {
            "cells": [[{"type": "text", "text": {"content": str(c)}}] for c in cells]
        }
    }


def _table_block(rows: list, has_header: bool = True) -> dict:
    """Cria um bloco tabela com linhas. rows[0] = cabeçalho se has_header=True."""
    return {
        "object": "block",
        "type": "table",
        "table": {
            "table_width": len(rows[0]) if rows else 1,
            "has_column_header": has_header,
            "has_row_header": False,
            "children": [_table_row(r) for r in rows]
        }
    }


# ── Setup: criação dos databases ─────────────────────────────────────────────

def setup_notion_databases(parent_page_id: str = NOTION_PAGE_DOCKING) -> dict:
    """
    Cria (ou recupera existentes) os três databases principais no Notion:
      - Experimentos: log de todas as análises rodadas
      - Compostos: cadastro de compostos com SMILES e metadados
      - Referências: literatura bibliográfica

    Retorna dict com os IDs: {"experiments": id, "compounds": id, "references": id}
    """
    # Verifica se já existem
    existing = _get_db_ids()
    if all(k in existing for k in ("experiments", "compounds", "references")):
        print(f"[Notion] Databases já configurados: {existing}")
        return existing

    print("[Notion] Criando databases no Notion...")
    db_ids = dict(existing)

    # ── Database de Experimentos ──────────────────────────────────────────
    if "experiments" not in db_ids:
        exp_db = _notion_request("POST", "databases", {
            "parent": {"type": "page_id", "page_id": parent_page_id},
            "icon": {"type": "emoji", "emoji": "🧬"},
            "title": [{"type": "text", "text": {"content": "📋 Experimentos — NQ-Docking"}}],
            "properties": {
                "Nome": {"title": {}},
                "Composto": {"rich_text": {}},
                "SMILES": {"rich_text": {}},
                "Tipo de Análise": {
                    "select": {
                        "options": [
                            {"name": "ADMET", "color": "blue"},
                            {"name": "Docking", "color": "green"},
                            {"name": "Toxicidade", "color": "red"},
                            {"name": "Ciclodextrina", "color": "yellow"},
                            {"name": "Dinâmica Molecular", "color": "purple"},
                            {"name": "Completo", "color": "default"},
                        ]
                    }
                },
                "Data": {"date": {}},
                "Session Tag": {"rich_text": {}},
                "MW (g/mol)": {"number": {"format": "number"}},
                "LogP": {"number": {"format": "number"}},
                "QED": {"number": {"format": "number"}},
                "Lipinski Pass": {
                    "select": {
                        "options": [
                            {"name": "Sim", "color": "green"},
                            {"name": "Não", "color": "red"},
                        ]
                    }
                },
                "GI Absorption": {
                    "select": {
                        "options": [
                            {"name": "High", "color": "green"},
                            {"name": "Low", "color": "orange"},
                        ]
                    }
                },
                "Score Docking (kcal/mol)": {"number": {"format": "number"}},
                "Alvo Proteico": {"rich_text": {}},
                "Status": {
                    "select": {
                        "options": [
                            {"name": "Concluído", "color": "green"},
                            {"name": "Em andamento", "color": "yellow"},
                            {"name": "Erro", "color": "red"},
                        ]
                    }
                },
                "Notas": {"rich_text": {}},
            }
        })
        db_ids["experiments"] = exp_db["id"]
        print(f"[Notion] ✅ Database Experimentos criado: {exp_db['id']}")

    # ── Database de Compostos ────────────────────────────────────────────
    if "compounds" not in db_ids:
        cmp_db = _notion_request("POST", "databases", {
            "parent": {"type": "page_id", "page_id": parent_page_id},
            "icon": {"type": "emoji", "emoji": "⚗️"},
            "title": [{"type": "text", "text": {"content": "⚗️ Compostos — NQ-Docking"}}],
            "properties": {
                "Nome": {"title": {}},
                "Abreviação": {"rich_text": {}},
                "SMILES": {"rich_text": {}},
                "Fórmula Molecular": {"rich_text": {}},
                "Peso Molecular": {"number": {"format": "number"}},
                "LogP": {"number": {"format": "number"}},
                "Classe": {
                    "select": {
                        "options": [
                            {"name": "Naftoquinona", "color": "orange"},
                            {"name": "Controle", "color": "blue"},
                            {"name": "Fármaco Aprovado", "color": "green"},
                            {"name": "Derivado", "color": "yellow"},
                        ]
                    }
                },
                "Ativo Antiparasitário": {
                    "select": {
                        "options": [
                            {"name": "Sim", "color": "green"},
                            {"name": "Não", "color": "red"},
                            {"name": "Pendente", "color": "gray"},
                        ]
                    }
                },
                "CAS": {"rich_text": {}},
                "PubChem CID": {"rich_text": {}},
                "Adicionado em": {"date": {}},
                "Origem": {"rich_text": {}},
            }
        })
        db_ids["compounds"] = cmp_db["id"]
        print(f"[Notion] ✅ Database Compostos criado: {cmp_db['id']}")

    # ── Database de Referências ──────────────────────────────────────────
    if "references" not in db_ids:
        ref_db = _notion_request("POST", "databases", {
            "parent": {"type": "page_id", "page_id": parent_page_id},
            "icon": {"type": "emoji", "emoji": "📚"},
            "title": [{"type": "text", "text": {"content": "📚 Referências — NQ-Docking"}}],
            "properties": {
                "Título": {"title": {}},
                "Autores": {"rich_text": {}},
                "Ano": {"number": {"format": "number"}},
                "Journal": {"rich_text": {}},
                "DOI": {"url": {}},
                "PubMed ID": {"rich_text": {}},
                "Tipo": {
                    "select": {
                        "options": [
                            {"name": "Artigo Original", "color": "blue"},
                            {"name": "Review", "color": "green"},
                            {"name": "Tese", "color": "yellow"},
                            {"name": "Livro", "color": "orange"},
                            {"name": "Preprint", "color": "gray"},
                        ]
                    }
                },
                "Relevância": {
                    "select": {
                        "options": [
                            {"name": "Alta", "color": "red"},
                            {"name": "Média", "color": "yellow"},
                            {"name": "Baixa", "color": "gray"},
                        ]
                    }
                },
                "Tema": {"multi_select": {
                    "options": [
                        {"name": "Docking Molecular", "color": "blue"},
                        {"name": "Naftoquinonas", "color": "orange"},
                        {"name": "Antiparasitário", "color": "green"},
                        {"name": "ADMET", "color": "purple"},
                        {"name": "Dinâmica Molecular", "color": "red"},
                        {"name": "Ciclodextrina", "color": "yellow"},
                    ]
                }},
                "Resumo": {"rich_text": {}},
                "Lido": {"checkbox": {}},
            }
        })
        db_ids["references"] = ref_db["id"]
        print(f"[Notion] ✅ Database Referências criado: {ref_db['id']}")

    _save_db_ids(db_ids)
    print(f"[Notion] ✅ Setup completo. Databases: {db_ids}")
    return db_ids


# ── Log automático de experimento ────────────────────────────────────────────

def log_experiment(
    session_tag: str,
    compound_name: str,
    analysis_type: str,          # "admet" | "docking" | "toxicity" | "cyclodextrin" | "md" | "complete"
    results: dict = None,
    notes: str = "",
    smiles: str = "",
) -> Optional[str]:
    """
    Registra uma análise no database de Experimentos do Notion.
    Retorna o ID da página criada, ou None em caso de erro.
    """
    db_ids = _get_db_ids()
    db_id = db_ids.get("experiments")
    if not db_id:
        print("[Notion] ⚠️  Database de Experimentos não configurado. Execute setup_notion_databases() primeiro.")
        return None

    results = results or {}
    tipo_map = {
        "admet": "ADMET",
        "docking": "Docking",
        "toxicity": "Toxicidade",
        "toxicidade": "Toxicidade",
        "cyclodextrin": "Ciclodextrina",
        "ciclodextrina": "Ciclodextrina",
        "md": "Dinâmica Molecular",
        "complete": "Completo",
    }
    tipo = tipo_map.get(analysis_type.lower(), "ADMET")

    # Extrai métricas principais do resultado
    mw   = float(results.get("MW", 0) or 0)
    logp = float(results.get("LogP", 0) or results.get("logP", 0) or 0)
    qed  = float(results.get("QED", 0) or results.get("qed", 0) or 0)
    lip  = results.get("Lipinski_pass", results.get("lipinski_pass", None))
    gi   = results.get("GI_absorption", results.get("gi_absorption", ""))
    score_docking = float(results.get("best_score", results.get("vina_score", 0)) or 0)
    alvo = results.get("target", results.get("receptor", ""))

    props = {
        "Nome": {"title": [{"text": {"content": f"{compound_name} — {tipo} — {session_tag}"}}]},
        "Composto": {"rich_text": [{"text": {"content": compound_name}}]},
        "Session Tag": {"rich_text": [{"text": {"content": session_tag}}]},
        "Tipo de Análise": {"select": {"name": tipo}},
        "Data": {"date": {"start": datetime.datetime.now().isoformat()}},
        "Status": {"select": {"name": "Concluído"}},
    }
    if smiles:
        props["SMILES"] = {"rich_text": [{"text": {"content": smiles[:2000]}}]}
    if mw:
        props["MW (g/mol)"] = {"number": mw}
    if logp:
        props["LogP"] = {"number": logp}
    if qed:
        props["QED"] = {"number": round(qed, 4)}
    if lip is not None:
        props["Lipinski Pass"] = {"select": {"name": "Sim" if lip else "Não"}}
    if gi:
        props["GI Absorption"] = {"select": {"name": gi if gi in ("High", "Low") else "High"}}
    if score_docking:
        props["Score Docking (kcal/mol)"] = {"number": score_docking}
    if alvo:
        props["Alvo Proteico"] = {"rich_text": [{"text": {"content": str(alvo)}}]}
    if notes:
        props["Notas"] = {"rich_text": [{"text": {"content": notes[:2000]}}]}

    try:
        page = _notion_request("POST", "pages", {
            "parent": {"database_id": db_id},
            "properties": props,
        })
        page_id = page["id"]
        print(f"[Notion] ✅ Experimento registrado: {page['url']}")
        return page_id
    except Exception as e:
        print(f"[Notion] ❌ Erro ao registrar experimento: {e}")
        return None


# ── Exportar relatório como página Notion ────────────────────────────────────

def export_report(
    session_tag: str,
    compound_name: str,
    admet_data: dict = None,
    tox_data: dict = None,
    docking_results: dict = None,
    text_en: str = "",
    text_pt: str = "",
    parent_page_id: str = NOTION_PAGE_DOCKING,
) -> Optional[str]:
    """
    Cria uma página de relatório completo no Notion com:
    - Resumo executivo (EN/PT)
    - Tabelas de métricas ADMET
    - Resultados de docking
    - Perfil de toxicidade
    Retorna URL da página criada.
    """
    admet_data    = admet_data    or {}
    tox_data      = tox_data      or {}
    docking_results = docking_results or {}

    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    title   = f"📄 Relatório {compound_name} — {session_tag}"

    # ── Constrói blocos do relatório ──────────────────────────────────────
    blocks = [
        _callout_block(
            f"Relatório gerado automaticamente pelo NQ-Docking Suite\n"
            f"Composto: {compound_name} | Sessão: {session_tag} | Data: {now_str}",
            "🧬"
        ),
        _heading_block("Executive Summary / Resumo Executivo", 2),
    ]

    if text_en:
        blocks.append(_heading_block("🇺🇸 English", 3))
        # Split long text into 2000-char paragraphs
        for chunk in [text_en[i:i+1900] for i in range(0, min(len(text_en), 9500), 1900)]:
            blocks.append(_text_block(chunk))

    if text_pt:
        blocks.append(_heading_block("🇧🇷 Português", 3))
        for chunk in [text_pt[i:i+1900] for i in range(0, min(len(text_pt), 9500), 1900)]:
            blocks.append(_text_block(chunk))

    blocks.append(_divider_block())

    # ── Tabela ADMET ──────────────────────────────────────────────────────
    if admet_data:
        blocks.append(_heading_block("🔬 Perfil ADMET", 2))
        admet_rows = [
            ["Propriedade", "Valor", "Classificação"],
            ["Peso Molecular (MW)", f"{admet_data.get('MW', 'N/A')} g/mol",
             "✅ OK" if float(admet_data.get("MW") or 0) < 500 else "⚠️ Alto"],
            ["LogP", str(admet_data.get("LogP", "N/A")),
             "✅ OK" if -0.4 <= float(admet_data.get("LogP") or 0) <= 5.6 else "⚠️ Fora"],
            ["HBD", str(admet_data.get("HBD", "N/A")),
             "✅" if int(admet_data.get("HBD") or 0) <= 5 else "⚠️"],
            ["HBA", str(admet_data.get("HBA", "N/A")),
             "✅" if int(admet_data.get("HBA") or 0) <= 10 else "⚠️"],
            ["TPSA", f"{admet_data.get('TPSA', 'N/A')} Å²",
             "✅ Oral" if float(admet_data.get("TPSA") or 0) < 140 else "⚠️ Alto"],
            ["QED", str(round(float(admet_data.get("QED") or 0), 3)), "—"],
            ["Lipinski", "✅ Passou" if admet_data.get("Lipinski_pass") else "❌ Falhou", "—"],
            ["GI Absorption", str(admet_data.get("GI_absorption", "N/A")), "—"],
            ["BBB Penetrant", "Sim" if admet_data.get("bbb_penetrant") else "Não", "—"],
            ["hERG Risk", str(admet_data.get("herg_risk", "N/A")), "—"],
            ["AMES Alert", "⚠️ Sim" if admet_data.get("ames_alert") else "✅ Não", "—"],
        ]
        blocks.append(_table_block(admet_rows))
        blocks.append(_divider_block())

    # ── Resultados de Docking ─────────────────────────────────────────────
    if docking_results:
        blocks.append(_heading_block("🎯 Resultados de Docking Molecular", 2))
        dock_rows = [["Alvo Proteico", "Melhor Score (kcal/mol)", "Modo"]]
        for target, data in docking_results.items():
            if isinstance(data, dict):
                score = data.get("best_score", data.get("score", "N/A"))
                mode  = data.get("mode", "1")
                dock_rows.append([target, str(score), str(mode)])
            else:
                dock_rows.append([target, str(data), "—"])
        if len(dock_rows) > 1:
            blocks.append(_table_block(dock_rows))
        blocks.append(_divider_block())

    # ── Toxicidade ───────────────────────────────────────────────────────
    if tox_data:
        blocks.append(_heading_block("☠️ Perfil de Toxicidade", 2))
        tox_rows = [["Parâmetro", "Valor"]]
        tox_keys = [
            ("ld50_oral_mg_kg", "LD50 Oral (mg/kg)"),
            ("ghs_class", "Classe GHS"),
            ("ames_positive", "AMES Mutagênico"),
            ("hepatotoxicity", "Hepatotoxicidade"),
            ("carcinogenicity", "Carcinogenicidade"),
            ("recommendation", "Recomendação"),
        ]
        for key, label in tox_keys:
            val = tox_data.get(key)
            if val is not None:
                tox_rows.append([label, str(val)])
        if len(tox_rows) > 1:
            blocks.append(_table_block(tox_rows))
        blocks.append(_divider_block())

    # ── Metadados da sessão ───────────────────────────────────────────────
    blocks.append(_heading_block("📁 Metadados da Sessão", 2))
    blocks.append(_code_block(
        f"session_tag: {session_tag}\n"
        f"compound: {compound_name}\n"
        f"generated: {now_str}\n"
        f"tool: NQ-Docking Suite v2 — IFRJ",
        "plain text"
    ))

    try:
        page = _notion_request("POST", "pages", {
            "parent": {"page_id": parent_page_id},
            "icon": {"type": "emoji", "emoji": "📄"},
            "properties": {
                "title": {"title": [{"text": {"content": title}}]}
            },
            "children": blocks[:100],  # Notion limit: 100 blocks per request
        })
        url = page.get("url", "")
        page_id = page["id"]
        print(f"[Notion] ✅ Relatório exportado: {url}")

        # Se há mais de 100 blocos, adiciona em lotes
        if len(blocks) > 100:
            remaining = blocks[100:]
            while remaining:
                batch = remaining[:100]
                remaining = remaining[100:]
                _notion_request("PATCH", f"blocks/{page_id}/children", {"children": batch})

        return url
    except Exception as e:
        print(f"[Notion] ❌ Erro ao exportar relatório: {e}")
        return None


# ── Sincronizar composto ─────────────────────────────────────────────────────

def sync_compound(
    compound_name: str,
    smiles: str = "",
    formula: str = "",
    mw: float = 0,
    logp: float = 0,
    compound_class: str = "Naftoquinona",
    antiparasitic: str = "Pendente",
    cas: str = "",
    pubchem_cid: str = "",
    source: str = "NQ-Docking Suite",
) -> Optional[str]:
    """
    Cria ou atualiza um composto no database de Compostos do Notion.
    Retorna o ID da página criada/atualizada.
    """
    db_ids = _get_db_ids()
    db_id = db_ids.get("compounds")
    if not db_id:
        print("[Notion] ⚠️  Database de Compostos não configurado. Execute setup_notion_databases() primeiro.")
        return None

    # Verifica se composto já existe
    try:
        existing = _notion_request("POST", f"databases/{db_id}/query", {
            "filter": {
                "property": "Nome",
                "title": {"equals": compound_name}
            }
        })
        if existing.get("results"):
            page_id = existing["results"][0]["id"]
            print(f"[Notion] Composto '{compound_name}' já existe, atualizando...")
        else:
            page_id = None
    except Exception:
        page_id = None

    props = {
        "Nome": {"title": [{"text": {"content": compound_name}}]},
        "Adicionado em": {"date": {"start": datetime.datetime.now().isoformat()}},
        "Classe": {"select": {"name": compound_class}},
        "Ativo Antiparasitário": {"select": {"name": antiparasitic}},
        "Origem": {"rich_text": [{"text": {"content": source}}]},
    }
    if smiles:
        props["SMILES"] = {"rich_text": [{"text": {"content": smiles[:2000]}}]}
    if formula:
        props["Fórmula Molecular"] = {"rich_text": [{"text": {"content": formula}}]}
    if mw:
        props["Peso Molecular"] = {"number": mw}
    if logp:
        props["LogP"] = {"number": logp}
    if cas:
        props["CAS"] = {"rich_text": [{"text": {"content": cas}}]}
    if pubchem_cid:
        props["PubChem CID"] = {"rich_text": [{"text": {"content": str(pubchem_cid)}}]}

    try:
        if page_id:
            page = _notion_request("PATCH", f"pages/{page_id}", {"properties": props})
        else:
            page = _notion_request("POST", "pages", {
                "parent": {"database_id": db_id},
                "properties": props,
            })
        print(f"[Notion] ✅ Composto sincronizado: {compound_name} → {page.get('url', '')}")
        return page["id"]
    except Exception as e:
        print(f"[Notion] ❌ Erro ao sincronizar composto: {e}")
        return None


# ── Importar compostos do Notion ─────────────────────────────────────────────

def import_compounds() -> list:
    """
    Importa compostos do database Notion para uso no projeto.
    Retorna lista de dicts com keys: name, smiles, formula, mw, logp, class.
    """
    db_ids = _get_db_ids()
    db_id = db_ids.get("compounds")
    if not db_id:
        print("[Notion] ⚠️  Database de Compostos não configurado.")
        return []

    try:
        result = _notion_request("POST", f"databases/{db_id}/query", {"page_size": 100})
        compounds = []
        for page in result.get("results", []):
            props = page.get("properties", {})

            def _get_title(prop_name):
                p = props.get(prop_name, {})
                t = p.get("title", [])
                return t[0]["text"]["content"] if t else ""

            def _get_text(prop_name):
                p = props.get(prop_name, {})
                rt = p.get("rich_text", [])
                return rt[0]["text"]["content"] if rt else ""

            def _get_number(prop_name):
                p = props.get(prop_name, {})
                return p.get("number") or 0

            def _get_select(prop_name):
                p = props.get(prop_name, {})
                s = p.get("select")
                return s["name"] if s else ""

            name = _get_title("Nome")
            if not name:
                continue
            compounds.append({
                "name":    name,
                "smiles":  _get_text("SMILES"),
                "formula": _get_text("Fórmula Molecular"),
                "mw":      _get_number("Peso Molecular"),
                "logp":    _get_number("LogP"),
                "class":   _get_select("Classe"),
                "abbr":    _get_text("Abreviação"),
                "cas":     _get_text("CAS"),
                "pubchem": _get_text("PubChem CID"),
                "notion_id": page["id"],
            })

        print(f"[Notion] ✅ {len(compounds)} compostos importados do Notion.")
        return compounds
    except Exception as e:
        print(f"[Notion] ❌ Erro ao importar compostos: {e}")
        return []


# ── Importar referências do Notion ───────────────────────────────────────────

def import_references() -> list:
    """
    Importa referências bibliográficas do database Notion.
    Retorna lista de dicts com metadados bibliográficos.
    """
    db_ids = _get_db_ids()
    db_id = db_ids.get("references")
    if not db_id:
        print("[Notion] ⚠️  Database de Referências não configurado.")
        return []

    try:
        result = _notion_request("POST", f"databases/{db_id}/query", {"page_size": 100})
        refs = []
        for page in result.get("results", []):
            props = page.get("properties", {})

            def _get_title(prop_name):
                p = props.get(prop_name, {})
                t = p.get("title", [])
                return t[0]["text"]["content"] if t else ""

            def _get_text(prop_name):
                p = props.get(prop_name, {})
                rt = p.get("rich_text", [])
                return rt[0]["text"]["content"] if rt else ""

            def _get_url(prop_name):
                p = props.get(prop_name, {})
                return p.get("url") or ""

            def _get_number(prop_name):
                p = props.get(prop_name, {})
                return p.get("number") or 0

            def _get_select(prop_name):
                p = props.get(prop_name, {})
                s = p.get("select")
                return s["name"] if s else ""

            def _get_multiselect(prop_name):
                p = props.get(prop_name, {})
                ms = p.get("multi_select", [])
                return [item["name"] for item in ms]

            title = _get_title("Título")
            if not title:
                continue
            refs.append({
                "title":   title,
                "authors": _get_text("Autores"),
                "year":    _get_number("Ano"),
                "journal": _get_text("Journal"),
                "doi":     _get_url("DOI"),
                "pmid":    _get_text("PubMed ID"),
                "type":    _get_select("Tipo"),
                "relevance": _get_select("Relevância"),
                "topics":  _get_multiselect("Tema"),
                "summary": _get_text("Resumo"),
                "notion_id": page["id"],
            })

        print(f"[Notion] ✅ {len(refs)} referências importadas do Notion.")
        return refs
    except Exception as e:
        print(f"[Notion] ❌ Erro ao importar referências: {e}")
        return []


# ── Sync em lote dos compostos do populate_ligand_db ────────────────────────

def sync_all_compounds_from_db() -> int:
    """
    Sincroniza todos os compostos cadastrados em populate_ligand_db.py para o Notion.
    Retorna o número de compostos sincronizados.
    """
    sys.path.insert(0, str(SCRIPT_DIR))
    try:
        from compound_db import list_compounds
        compounds = list_compounds()
    except Exception as e:
        print(f"[Notion] ❌ Erro ao carregar compostos locais: {e}")
        return 0

    count = 0
    for c in compounds:
        name   = c.get("name", "")
        smiles = c.get("smiles", "")
        if not name:
            continue
        ok = sync_compound(
            compound_name   = name,
            smiles          = smiles,
            formula         = c.get("formula", ""),
            mw              = float(c.get("mw") or 0),
            logp            = float(c.get("logp") or 0),
            compound_class  = "Naftoquinona" if "naphth" in name.lower() or "nq" in name.lower() else "Controle",
            antiparasitic   = "Pendente",
            source          = "populate_ligand_db.py",
        )
        if ok:
            count += 1

    print(f"[Notion] ✅ {count}/{len(compounds)} compostos sincronizados.")
    return count


# ── Status / health check ────────────────────────────────────────────────────

def get_status() -> dict:
    """Retorna status da integração Notion (token, databases configurados, etc.)."""
    token = _get_token()
    db_ids = _get_db_ids()
    status = {
        "token_configured": bool(token),
        "token_preview": f"{token[:12]}..." if token else None,
        "databases": db_ids,
        "databases_ready": all(k in db_ids for k in ("experiments", "compounds", "references")),
        "config_file": str(CONFIG_FILE),
        "config_exists": CONFIG_FILE.exists(),
    }
    if token:
        try:
            me = _notion_request("GET", "users/me")
            status["workspace"] = me.get("bot", {}).get("workspace_name", "")
            status["bot_name"]  = me.get("name", "")
            status["connected"] = True
        except Exception as e:
            status["connected"] = False
            status["error"] = str(e)
    else:
        status["connected"] = False

    return status


# ── CLI standalone ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="NQ-Docking ↔ Notion Integration")
    parser.add_argument("command", choices=["setup", "status", "sync_compounds", "import_compounds"],
                        help="Comando a executar")
    args = parser.parse_args()

    if args.command == "setup":
        setup_notion_databases()
    elif args.command == "status":
        s = get_status()
        print(json.dumps(s, indent=2, ensure_ascii=False))
    elif args.command == "sync_compounds":
        n = sync_all_compounds_from_db()
        print(f"Sincronizados: {n}")
    elif args.command == "import_compounds":
        comps = import_compounds()
        for c in comps:
            print(f"  {c['name']} | {c['smiles'][:40]}...")
