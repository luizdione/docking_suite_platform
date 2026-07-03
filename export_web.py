#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""export_web.py - empacota um BUNDLE ESTATICO publicavel (NQ-Docking Viewer).

Gera um diretorio 100% estatico (HTML/JS/JSON) que reproduz o dashboard em modo
**somente leitura** + uma calculadora ADMET que roda no navegador, SEM nenhuma
dependencia de Python/GPU/rede no servidor. E o que se hospeda na Hostinger em
`docking.ifrj-crj-geneticamolecular.online` (par do bonder/epigenoma).

Ideia (igual ao bonder): a analise pesada roda LOCAL (dashboard em localhost:5001,
docking com Vina, MD com OpenMM/GPU). Este exportador captura os RESULTADOS ja
gerados pelas rotas de leitura do dashboard (via test_client), **sanitiza os
caminhos absolutos locais** (nunca publica C:\\... ou E:\\...), copia so imagens
leves (nunca .dcd/.chk/.pdb pesados) e monta o pacote para envio por FTP/hPanel.

Estrutura do bundle:
  index.html                viewer estatico (web/viewer/index.html)
  lib/*.js                  marked (vendorizado; sem CDN)
  .htaccess                 MIME + no-indexes (+ Basic Auth comentado)
  README_DEPLOY.md          passo a passo do deploy na Hostinger
  config.json               versao, composto, alvos, seletividade, subprojeto
  data/compounds.json       biblioteca de compostos (sanitizada)
  data/admet.json           perfil ADMET pre-computado
  data/binding_mode.json    modos subversivo vs inibidor (sanitizado)
  data/md_index.json        tags de MD + figuras RMSD disponiveis
  data/md/<tag>.json        resumo por corrida de MD
  data/subprojects.json     guarda-chuva (programa + subprojetos)
  figures/*.png|svg         imagens leves (ADMET, estrutura, RMSD)
  assets/pipeline_SUB01.svg SVG do pipeline do subprojeto ativo
  docs/index.json + *.md    documentos da aba Ajuda

Uso:
  python export_web.py --out dist_web
  python export_web.py --out dist_web --compound 2-SCN-NQ --max-img-mb 3
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
VIEWER_DIR = SCRIPT_DIR / "web" / "viewer"
BASE_DIR = Path(os.environ.get("NQ_DATA_DIR", r"E:\data\docking"))

# chaves com caminhos/flags locais que NAO devem ir para o site publico
DROP_KEYS = {
    "_dir", "_has_mol2", "_has_pdbqt", "path", "docking_dir", "figures_dir",
    "logs_dir", "results_file", "texts_dir", "kit_dir",
}
IMG_EXT = {".png", ".svg", ".jpg", ".jpeg", ".webp", ".gif"}
_ABS_PATH = re.compile(r"^[A-Za-z]:[\\/]|^/(home|Users|mnt|data)/")

# Painel de alvos (PDB/organismo/categoria/cofatores) — identidade visual do dashboard.
# category: parasite (T. cruzi) | leishmania | human — controla a cor do badge PDB.
TARGET_INFO = [
    {"code": "TcTR", "pdb": "1BZL", "organism": "Trypanosoma cruzi", "name": "Trypanothione Reductase", "category": "parasite", "cofactors": ["FAD"]},
    {"code": "TcCruzain", "pdb": "3KKU", "organism": "Trypanosoma cruzi", "name": "Cruzain", "category": "parasite", "cofactors": []},
    {"code": "TcFeSOD", "pdb": "4DVH", "organism": "Trypanosoma cruzi", "name": "Fe-Superoxide Dismutase", "category": "parasite", "cofactors": ["FE", "FE2"]},
    {"code": "LiTR", "pdb": "2JK6", "organism": "Leishmania infantum (proxy L. amazonensis)", "name": "Trypanothione Reductase", "category": "leishmania", "cofactors": ["FAD"]},
    {"code": "LmPTR1", "pdb": "2QHX", "organism": "Leishmania major", "name": "Pteridine Reductase 1", "category": "leishmania", "cofactors": ["NAP", "NDP", "FAD"]},
    {"code": "LmDHODH", "pdb": "3C61", "organism": "Leishmania donovani", "name": "Dihydroorotate Dehydrogenase", "category": "leishmania", "cofactors": ["FMN", "FAD"]},
    {"code": "hGR", "pdb": "3DJJ", "organism": "Homo sapiens", "name": "Glutathione Reductase", "category": "human", "cofactors": ["FAD"]},
    {"code": "hCatL", "pdb": "2XU3", "organism": "Homo sapiens", "name": "Cathepsin L", "category": "human", "cofactors": []},
    {"code": "hDHODH", "pdb": "1D3H", "organism": "Homo sapiens", "name": "Dihydroorotate Dehydrogenase", "category": "human", "cofactors": ["FMN"]},
    {"code": "hMnSOD", "pdb": "1LUV", "organism": "Homo sapiens", "name": "Mn-Superoxide Dismutase", "category": "human", "cofactors": ["MN"]},
]


def _looks_like_local_path(s: str) -> bool:
    return bool(_ABS_PATH.match(s)) or ("\\" in s and ":" in s)


def sanitize(obj):
    """Remove chaves de path e strings que sejam caminhos absolutos locais."""
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items() if k not in DROP_KEYS}
    if isinstance(obj, list):
        return [sanitize(x) for x in obj]
    if isinstance(obj, str) and _looks_like_local_path(obj):
        return None
    return obj


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def _get_json(client, url):
    r = client.get(url)
    if r.status_code != 200:
        return None
    try:
        return json.loads(r.get_data(as_text=True))
    except Exception:
        return None


def _copy_image(src: Path, figures_out: Path, max_bytes: int) -> str | None:
    """Copia uma imagem leve para figures/ (nome unico). Retorna o rel path ou None."""
    if not src.exists() or src.suffix.lower() not in IMG_EXT:
        return None
    if src.stat().st_size > max_bytes:
        return None
    figures_out.mkdir(parents=True, exist_ok=True)
    dest = figures_out / src.name
    if dest.exists() and dest.stat().st_size == src.stat().st_size:
        return f"figures/{dest.name}"
    # evita colisao de nomes iguais com conteudos diferentes
    if dest.exists():
        dest = figures_out / f"{src.stem}_{abs(hash(str(src))) % 9999}{src.suffix}"
    shutil.copyfile(src, dest)
    return f"figures/{dest.name}"


def export_bundle(out_dir: str, compound: str = "2-SCN-NQ", max_img_mb: float = 3.0) -> dict:
    out = Path(out_dir)
    max_bytes = int(max_img_mb * 1024 * 1024)

    # 1) viewer estatico + libs + assets de deploy (copia a arvore inteira)
    if not VIEWER_DIR.exists():
        raise FileNotFoundError(f"viewer nao encontrado em {VIEWER_DIR}")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(VIEWER_DIR, out, dirs_exist_ok=True)

    figures_out = out / "figures"
    summary = {"out_dir": str(out), "compound": compound, "figures": 0, "md_tags": 0, "docs": 0}

    # 2) sobe o dashboard em test_client (rotas de leitura; SEM GPU)
    sys.path.insert(0, str(SCRIPT_DIR))
    import dashboard_v2 as D
    client = D.app.test_client()

    # --- compostos ---
    comp = _get_json(client, "/compounds") or {}
    _write_json(out / "data" / "compounds.json", sanitize(comp))

    # --- ADMET (arquivo pre-computado tem prioridade; senao recalcula on-the-fly) ---
    admet = None
    admet_file = BASE_DIR / "results" / compound / "admet_profile.json"
    if admet_file.exists():
        try:
            raw = json.loads(admet_file.read_text(encoding="utf-8"))
            admet = raw if "admet" in raw else {"admet": raw}
        except Exception:
            admet = None
    if admet is None:
        admet = _get_json(client, "/admet_profile") or {"admet": {}}
    _write_json(out / "data" / "admet.json", sanitize(admet))

    # --- binding mode (subversivo vs inibidor) ---
    binding_raw = _get_json(client, "/binding_mode_status") or {}
    _write_json(out / "data" / "binding_mode.json", sanitize(binding_raw))

    # --- MD: tags + resumo por tag + figuras RMSD ---
    tags = (_get_json(client, "/md_rmsd_tags") or {}).get("tags", [])
    md_figs = []
    for tag in tags:
        s = _get_json(client, f"/md_rmsd_summary/{tag}")
        if s is not None:
            _write_json(out / "data" / "md" / f"{tag}.json", sanitize(s))
    # figuras RMSD: procura *rmsd*.png em figures/ (E:) e docking_project/figures/
    search_roots = [BASE_DIR / "figures", SCRIPT_DIR / "docking_project" / "figures"]
    seen = set()
    for root in search_roots:
        if not root.exists():
            continue
        for img in root.rglob("*.png"):
            low = img.name.lower()
            if "rmsd" not in low and "rmsf" not in low:
                continue
            if img.name in seen:
                continue
            rel = _copy_image(img, figures_out, max_bytes)
            if rel:
                seen.add(img.name)
                # casa a figura a uma tag por substring, preferindo a tag mais especifica (longa)
                tag_match = next((t for t in sorted(tags, key=len, reverse=True)
                                  if t.lower() in low), None)
                md_figs.append({"name": img.name, "rel": rel, "tag": tag_match})
    _write_json(out / "data" / "md_index.json", {"tags": tags, "figures": md_figs})
    summary["md_tags"] = len(tags)

    # --- figuras do composto ativo (ADMET, estrutura) ---
    figs = (_get_json(client, "/figures") or {}).get("figures", [])
    n_img = len(md_figs)
    for f in figs:
        p = Path(f.get("path", ""))
        if _copy_image(p, figures_out, max_bytes):
            n_img += 1
    summary["figures"] = n_img

    # --- subprojetos + pipeline SVG ---
    sp = _get_json(client, "/api/subprojects") or {}
    _write_json(out / "data" / "subprojects.json", sanitize(sp))
    active_sid = None
    for s in sp.get("subprojetos", []):
        active_sid = s.get("id")
        break
    if active_sid:
        r = client.get(f"/api/subproject/{active_sid}/pipeline.svg")
        if r.status_code == 200:
            (out / "assets").mkdir(parents=True, exist_ok=True)
            (out / "assets" / f"pipeline_{active_sid}.svg").write_bytes(r.get_data())
            # o viewer procura assets/pipeline_SUB01.svg
            if active_sid != "SUB01":
                shutil.copyfile(out / "assets" / f"pipeline_{active_sid}.svg",
                                out / "assets" / "pipeline_SUB01.svg")

    # --- docs (.md da aba Ajuda) ---
    docs_index = []
    help_list = (_get_json(client, "/help_list") or {}).get("docs", [])
    (out / "docs").mkdir(parents=True, exist_ok=True)
    for d in help_list:
        fname = d.get("file")
        if not fname:
            continue
        key = fname[:-3] if fname.lower().endswith(".md") else fname
        doc = _get_json(client, f"/help_doc/{key}")
        raw = (doc or {}).get("raw", "")
        if not raw:
            continue
        (out / "docs" / fname).write_text(raw, encoding="utf-8")
        docs_index.append({"id": key, "title": d.get("title", key), "file": fname})
    _write_json(out / "docs" / "index.json", {"docs": docs_index})
    summary["docs"] = len(docs_index)

    # --- config.json (contexto para o viewer) ---
    version = "0.0.0"
    vfile = SCRIPT_DIR / "VERSION"
    if vfile.exists():
        version = vfile.read_text(encoding="utf-8").strip()
    subproject = sp.get("subprojetos", [{}])[0] if sp.get("subprojetos") else {}
    # seletividade: distancia lig-FAD por alvo (modo subversivo), agregando n replicas
    selectivity = _build_selectivity(binding_raw)
    # figura de estrutura 2D (hero da Visao geral), se copiada para figures/
    structure_svg = None
    if figures_out.exists():
        for cand in sorted(figures_out.glob("*structure*.svg")) + sorted(figures_out.glob("*struct*.svg")):
            structure_svg = f"figures/{cand.name}"
            break
    config = {
        "title": "Docking Suite IFRJ · naftoquinonas antiparasitárias (docking + MD + ADMET)",
        "version": version,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "compound": compound,
        "paper_title": subproject.get("titulo_paper"),
        "structure_svg": structure_svg,
        "targets": TARGET_INFO,
        "selectivity": selectivity,
        "subproject": {
            "id": subproject.get("id"),
            "titulo_paper": subproject.get("titulo_paper"),
            "periodico_alvo": subproject.get("periodico_alvo"),
            "fase_atual": subproject.get("fase_atual"),
            "status": subproject.get("status"),
        },
    }
    _write_json(out / "config.json", config)

    return summary


def _build_selectivity(binding_raw):
    """Distancia lig-FAD (modo subversivo) por alvo, agregando as replicas n=1..3.

    Le os binding_mode_<alvo>_subversive[/_r2/_r3].json do kit_dir (analise 3D do
    projeto do paper, referenciado pela rota /binding_mode_status). Reproduz o achado
    central: TR parasitaria (baixa dist. FAD) vs GR humana (alta dist.). Silencioso se
    o kit_dir nao existir nesta maquina (selectivity fica []; o viewer lida com isso)."""
    kit = (binding_raw or {}).get("kit_dir")
    if not kit:
        return []
    kitp = Path(kit)
    if not kitp.exists():
        return []
    out = []
    for tgt, is_parasite in (("TcTR", True), ("LiTR", True), ("hGR", False)):
        vals = []
        for suffix in ("_subversive", "_subversive_50ns", "_subversive_r2", "_subversive_r3"):
            f = kitp / f"binding_mode_{tgt}{suffix}.json"
            if not f.exists():
                continue
            try:
                d = json.loads(f.read_text(encoding="utf-8"))
                m = (d.get("lig_fad_A") or {}).get("mean")
                if m is not None:
                    vals.append(float(m))
            except Exception:
                continue
        if vals:
            out.append({
                "target": tgt,
                "fad_dist": round(statistics.mean(vals), 2),
                "fad_sd": round(statistics.pstdev(vals), 2) if len(vals) > 1 else 0.0,
                "n": len(vals),
                "parasite": is_parasite,
            })
    return out


def main():
    ap = argparse.ArgumentParser(description="Empacota o NQ-Docking Viewer estatico (Hostinger).")
    ap.add_argument("--out", default="dist_web", help="diretorio de saida do bundle")
    ap.add_argument("--compound", default="2-SCN-NQ", help="composto base do bundle")
    ap.add_argument("--max-img-mb", type=float, default=3.0, help="tamanho maximo por imagem (MB)")
    args = ap.parse_args()
    s = export_bundle(args.out, compound=args.compound, max_img_mb=args.max_img_mb)
    print("[OK] bundle gerado:")
    for k, v in s.items():
        print(f"     {k}: {v}")
    print(f"\nProximo passo: envie o conteudo de {s['out_dir']}/ para o subdominio "
          "docking.ifrj-crj-geneticamolecular.online (ver README_DEPLOY.md).")


if __name__ == "__main__":
    main()
