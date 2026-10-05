

#!/usr/bin/env python3
"""
Dashboard v2 - Docking + ADMET + MD + Ciclodextrina + Toxicidade + Produção de Texto
Run: python dashboard_v2.py   |   Open: http://localhost:5001
FIX: usa substituicao direta de string em vez de Jinja2 para evitar conflito
     com template literals do JavaScript (${...})
"""
import os, sys, json, subprocess, io
from pathlib import Path

try:
    from flask import Flask, jsonify, request, Response
except ImportError:
    import _deps
    _deps.ensure(["flask"])  # raise se nao opt-in; instala se DOCKING_AUTO_PIP=1
    from flask import Flask, jsonify, request, Response

app        = Flask(__name__)
# Diretório base configurável via variável de ambiente NQ_DATA_DIR
# Padrão novo: E:\data\docking  (fallback: docking_project para compatibilidade)
import os as _os
BASE_DIR   = Path(_os.environ.get("NQ_DATA_DIR", r"E:\data\docking"))
SCRIPT_DIR = Path(__file__).resolve().parent


def _ensure_script_path():
    """Garante SCRIPT_DIR no sys.path SEM duplicar entradas a cada request.

    Antes, cada rota chamava sys.path.insert(0, ...) por requisicao, fazendo
    o sys.path crescer indefinidamente num servidor de longa duracao.
    """
    p = str(SCRIPT_DIR)
    if p not in sys.path:
        sys.path.insert(0, p)


# Insere uma vez no import (as rotas ainda chamam _ensure_script_path, agora idempotente)
_ensure_script_path()


def _clean_compound(value, default="2-SCN-NQ"):
    """Normaliza o nome do composto vindo do frontend.

    Rede de defesa contra objetos JS serializados ('[object Object]'),
    valores nao-string e vazios — que geravam pastas espurias '_object_Object_'
    via compound_db._safe_name. Retorna o composto default do projeto nesses casos.
    """
    if not isinstance(value, str):
        return default
    s = value.strip()
    if not s or s.lower() in ("[object object]", "undefined", "null"):
        return default
    return s

# ── HTML ──────────────────────────────────────────────────────────
# Usa __TARGETS__ como placeholder — SEM Jinja2 — para evitar que o
# motor de template tente processar os ${...} do JavaScript.
# ──────────────────────────────────────────────────────────────────
# ─── UI ───────────────────────────────────────────────────────────
# A interface vive em static/app.html (extraida de uma string de ~275 mil
# caracteres que ficava aqui). E servida CRUA, sem motor de template: Jinja
# colidiria com os template literals `${...}` do JS embutido. A injecao de
# dados continua sendo substituicao de placeholder, como antes.
_STATIC_DIR = Path(__file__).resolve().parent / "static"
_APP_HTML = _STATIC_DIR / "app.html"
_HTML_CACHE = {"mtime": None, "text": None}


def _load_html():
    """Le a UI do disco, com cache invalidado por mtime.

    Recarregar em desenvolvimento sem reiniciar o servidor e o ganho pratico
    imediato da extracao.
    """
    try:
        mtime = _APP_HTML.stat().st_mtime_ns
    except OSError as exc:
        raise RuntimeError(
            f"UI nao encontrada em {_APP_HTML}. O arquivo static/app.html faz parte "
            "do repositorio; se sumiu, restaure-o do git."
        ) from exc
    if _HTML_CACHE["mtime"] != mtime:
        _HTML_CACHE["text"] = _APP_HTML.read_text(encoding="utf-8")
        _HTML_CACHE["mtime"] = mtime
    return _HTML_CACHE["text"]


class _HTMLProxy(str):
    """Compatibilidade: o codigo existente usa `HTML.replace(...)`.

    Mantem `HTML` como se fosse a string de antes, mas relendo do disco a cada
    acesso. Assim nenhum chamador precisa mudar.
    """

    def __new__(cls):
        return super().__new__(cls, "")

    def replace(self, *args, **kwargs):
        return _load_html().replace(*args, **kwargs)

    def __str__(self):
        return _load_html()

    def __len__(self):
        return len(_load_html())

    def __contains__(self, item):
        return item in _load_html()


HTML = _HTMLProxy()



# ─── TARGET DB LOADER ─────────────────────────────────────────────
def _targets_json(campaign_id=None):
    """Alvos servidos ao front, vindos do registry da campanha.

    Havia aqui um terceiro dicionario de alvos hard-coded como fallback, e ele
    divergia da fonte real: TcCruzain apontava para 3GRS (que e glutationa
    redutase HUMANA), TcFeSOD para 1E7W e hDHODH para 4DVH (que e o Fe-SOD do
    parasito). Como o fallback so aparecia quando o import falhava, o erro
    passava despercebido e a UI exibia o PDB errado. A fonte agora e uma so.
    """
    targets_safe = {}
    try:
        _ensure_script_path()
        import registry

        for target in registry.list_targets(campaign_id or _default_campaign()):
            struct = target.get("structure") or {}
            chem = target.get("chemistry") or {}
            targets_safe[target["id"]] = {
                "pdb":            struct.get("pdb_id") or "N/A",
                "organism":       target.get("organism", ""),
                "target":         target.get("name", ""),
                "category":       target.get("category", "parasite"),
                "keep_cofactors": chem.get("keep_cofactors", []),
                "notes":          target.get("notes", ""),
                "structure_source": struct.get("source", "pdb"),
            }
    except Exception as exc:
        # Sem registry, ainda ha o fallback embutido de docking_v2 -- uma fonte
        # a menos que antes, e a que o teste de paridade vigia.
        import sys as _sys
        print(f"[AVISO] registry indisponivel ({exc}); caindo em docking_v2.TARGET_DB",
              file=_sys.stderr, flush=True)
        try:
            from docking_v2 import TARGET_DB
            for tid, info in TARGET_DB.items():
                targets_safe[tid] = {
                    "pdb":            info.get("pdb", "N/A"),
                    "organism":       info.get("organism", ""),
                    "target":         info.get("target", ""),
                    "category":       info.get("category", "parasite"),
                    "keep_cofactors": info.get("keep_cofactors", []),
                    "notes":          info.get("notes", ""),
                    "structure_source": "pdb",
                }
        except Exception as exc2:
            import traceback as _tb
            print(f"[ERRO] nenhuma fonte de alvos disponivel: {exc2}", file=_sys.stderr)
            _tb.print_exc()
    return targets_safe


def _default_campaign():
    """Campanha corrente. Default do programa, com escape por variavel de ambiente."""
    import os as _os

    env = _os.environ.get("DOCKING_CAMPAIGN")
    if env:
        return env
    try:
        import registry
        return registry.load_program().get("default_campaign", "SUB01")
    except Exception:
        return "SUB01"


# ─── ROUTES ───────────────────────────────────────────────────────

@app.route('/')
def index():
    # ?campaign= torna o dashboard multi-campanha (Fase 5.4); default = campanha do programa.
    campaign = request.args.get('campaign') or _default_campaign()
    targets_safe = _targets_json(campaign)
    # Bootstrap data-driven (Fase 5.2/5.4): campanha, ui, lista de campanhas p/ o seletor.
    try:
        _ensure_script_path()
        import web_api
        boot = web_api.bootstrap_for(campaign)
    except Exception as _e:  # nunca deixa o index cair por causa do bootstrap
        boot = {"campaign_id": campaign, "campaigns": [], "ui": {}, "error": str(_e)}
    # SUBSTITUICAO DIRETA — SEM Jinja2 — evita conflito com ${...} do JS
    page = (HTML.replace('__TARGETS__', json.dumps(targets_safe))
                .replace('__BOOTSTRAP__', json.dumps(boot)))
    return Response(page, mimetype='text/html')


@app.route('/run', methods=['POST'])
def run_pipeline():
    data            = request.get_json(silent=True) or {}
    targets         = data.get('targets', [])
    exh             = data.get('exhaustiveness', 32)
    compound_name   = _clean_compound(data.get('compound_name'))
    compound_smiles = data.get('compound_smiles', None)
    session_tag     = data.get('session_tag', None)

    log_buf = io.StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = log_buf
    success = False; actual_session = session_tag
    try:
        _ensure_script_path()
        from importlib import reload
        import docking_v2; reload(docking_v2)
        result = docking_v2.run_full_pipeline(
            targets,
            exhaustiveness=exh,
            compound_name=compound_name,
            compound_smiles=compound_smiles,
            session_tag=session_tag,
        )
        # run_full_pipeline now returns (bool, session_tag)
        if isinstance(result, tuple):
            success, actual_session = result
        else:
            success = bool(result)
    except Exception as e:
        import traceback
        log_buf.write(f"\n[ERRO] {e}\n{traceback.format_exc()}")
    finally:
        sys.stdout, sys.stderr = old_out, old_err

    # Load results from compound-specific path
    from compound_db import _safe_name
    results = {}
    rf = BASE_DIR / "results" / _safe_name(compound_name) / "results.json"
    if rf.exists():
        try:
            with open(rf) as f: results = json.load(f).get("results", {})
        except Exception as e:
            app.logger.warning("run_pipeline: falha ao ler %s: %s", rf, e)

    return jsonify({
        "success":     success,
        "log":         log_buf.getvalue(),
        "results":     results,
        "session_tag": actual_session,
        "compound":    compound_name,
    })


@app.route('/run_step', methods=['POST'])
def run_step():
    data = request.get_json(silent=True) or {}
    targets = data.get('targets', [])
    step = data.get('step', 'qc')
    log_buf = io.StringIO()
    old_out = sys.stdout; sys.stdout = log_buf
    try:
        if step == 'qc':
            from docking_v2 import TARGET_DB, download_pdb, analyze_pdb
            pd = BASE_DIR / "proteins"; pd.mkdir(parents=True, exist_ok=True)
            for tid in targets:
                if tid in TARGET_DB:
                    info = TARGET_DB[tid]
                    pdb_path, msg = download_pdb(info["pdb"], pd)
                    log_buf.write(f"\n{tid} ({info['pdb']}): {msg}\n")
                    if pdb_path:
                        qc = analyze_pdb(pdb_path, info)
                        for k in ['cofactors_found','metals_found','ligands_found','status']:
                            log_buf.write(f"  {k}: {qc.get(k)}\n")
        elif step == 'analysis':
            compound_name = _clean_compound(data.get('compound_name'))
            session_tag   = data.get('session_tag') or None
            from compound_db import get_session_dirs, make_session_tag, _safe_name
            if session_tag is None:
                session_tag = make_session_tag()
            dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)
            rf = dirs["results"] / "results.json"
            if rf.exists():
                from analysis import run_analysis
                for lang in ['en','pt']:
                    run_analysis(str(rf), str(dirs["figures"]), lang)
                log_buf.write(f"\n[OK] Analise completa — sessao: {session_tag}\n")
            else:
                log_buf.write("\n[ERRO] Sem results.json. Rode o docking primeiro.\n")
    except Exception as e:
        log_buf.write(f"\n[ERRO] {e}\n")
    finally:
        sys.stdout = old_out
    return jsonify({"log": log_buf.getvalue()})


@app.route('/run_admet', methods=['POST'])
def run_admet_route():
    data            = request.get_json(silent=True) or {}
    smiles          = data.get('smiles') or "O=C1C=C(SC#N)C(=O)c2ccccc21"
    lang            = data.get('lang', 'en')
    compound_name   = _clean_compound(data.get('compound_name'))
    session_tag     = data.get('session_tag') or None

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        session_tag = make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)

    log_buf = io.StringIO(); old_out = sys.stdout; sys.stdout = log_buf
    admet_flat = None; text_out = ""
    try:
        from importlib import reload
        import admet_module; reload(admet_module)
        figs_dir = dirs["figures"]
        # Pass compound_name so figure titles show the selected compound
        admet_full = admet_module.run_admet(smiles, str(figs_dir), lang,
                                            compound_name=compound_name,
                                            compound_abbr=compound_name)
        if lang != 'pt':
            admet_module.run_admet(smiles, str(figs_dir), 'pt',
                                   compound_name=compound_name,
                                   compound_abbr=compound_name)
        text_out = admet_module.generate_admet_text(admet_full, lang=lang)
        # Save text to session texts dir
        txt_file = dirs["texts"] / f"admet_{lang}.txt"
        txt_file.write_text(text_out, encoding='utf-8')
        if admet_full and "error" not in admet_full:
            pc = admet_full.get("physicochemical", {})
            ab = admet_full.get("absorption", {})
            dl = admet_full.get("drug_likeness", {})
            ex = admet_full.get("excretion", {})
            di = admet_full.get("distribution", {})
            me = admet_full.get("metabolism", {})
            to = admet_full.get("toxicity", {})
            admet_flat = {
                "MW": pc.get("MW",215.23), "LogP": pc.get("logP",2.41),
                "HBD": pc.get("HBD",0), "HBA": pc.get("HBA",3),
                "TPSA": pc.get("TPSA",67.72), "RotBonds": pc.get("rotatable_bonds",2),
                "QED": dl.get("qed",0.61),
                "GI_absorption": ab.get("oral_class","High"),
                "Caco2_Papp_nm_s": round(ab.get("caco2_papp",0)*1e7, 2),
                "Caco2_class": ab.get("caco2_class","High"),
                "logS_ESOL": ex.get("logS_esol",-3.3),
                "Solubility_class": ex.get("solubility_class","Soluble"),
                "Lipinski_pass": dl.get("lipinski",{}).get("pass",True),
                "Veber_pass": dl.get("veber",{}).get("pass",True),
                "bioavailability_score": ab.get("bioavailability_score",0.55),
                "hia_percent": ab.get("hia_percent",80),
                "bbb_penetrant": di.get("bbb_penetrant",False),
                "ppb_class": di.get("ppb_class","High (>90%)"),
                "t_half_hours": me.get("t_half_hours",4.0),
                "sas": dl.get("sas",2.5),
                "cyp_alerts": me.get("cyp_alerts",[]),
                "herg_risk": to.get("herg_risk","Low"),
                "ames_alert": to.get("ames_alert",False),
                # NOTA: LD50/Toxicity_class e CD_Ka_M1 NAO sao previstos aqui —
                # dominios proprios (toxicity_module / cyclodex_module). Nao gravar
                # constantes fabricadas no admet_profile.json (integridade de dados).
            }
            results_dir = dirs["results"]; results_dir.mkdir(parents=True, exist_ok=True)
            with open(results_dir / "admet_profile.json", "w") as fh:
                json.dump(admet_flat, fh, indent=2)
    except Exception as e:
        import traceback
        log_buf.write(f"\n[ERRO ADMET] {e}\n{traceback.format_exc()}")
    finally:
        sys.stdout = old_out
    return jsonify({"log": log_buf.getvalue(), "admet": admet_flat, "text": text_out,
                    "session_tag": session_tag, "figures_dir": str(dirs["figures"])})


@app.route('/admet_profile')
def admet_profile_route():
    compound = request.args.get('compound', '')
    # Search compound-specific results dir first
    if compound:
        from compound_db import _safe_name
        jp = BASE_DIR / "results" / _safe_name(compound) / "admet_profile.json"
        if jp.exists():
            with open(jp) as f: return jsonify({"admet": json.load(f)})
    # Fallback: scan all compound results dirs for any admet_profile.json
    if (BASE_DIR/"results").exists():
        for subdir in sorted((BASE_DIR/"results").iterdir(), reverse=True):
            jp = subdir / "admet_profile.json"
            if jp.exists():
                with open(jp) as f: return jsonify({"admet": json.load(f)})
    return jsonify({"admet": None})


@app.route('/figures')
def list_figures():
    """List figures — prefers compound+session path, falls back to legacy flat path."""
    compound    = request.args.get('compound', '')
    session_tag = request.args.get('session', '')
    figs = []

    if compound and session_tag:
        _ensure_script_path()
        from compound_db import get_session_dirs
        dirs = get_session_dirs(BASE_DIR, session_tag, compound)
        fd = dirs["figures"]
    else:
        # Legacy / fallback: scan all compound session dirs and return latest
        fd = None
        best_time = ""
        for subdir in sorted((BASE_DIR / "figures").iterdir()) if (BASE_DIR/"figures").exists() else []:
            if subdir.is_dir():
                for sess_dir in sorted(subdir.iterdir(), reverse=True):
                    if sess_dir.is_dir() and any(sess_dir.glob("*.png")):
                        if sess_dir.name > best_time:
                            best_time = sess_dir.name
                            fd = sess_dir
                        break

    if fd and fd.exists():
        for f in sorted(fd.iterdir()):
            if f.suffix in ['.png','.svg','.tiff']:
                figs.append({"name": f.name, "size": f.stat().st_size,
                             "path": str(f)})
    return jsonify({"figures": figs})


@app.route('/fig_img/<path:fname>')
def serve_fig(fname):
    """Serve a figure image — searches all compound/session subdirectories."""
    figs_base = BASE_DIR / "figures"
    base_resolved = str(figs_base.resolve())
    # Direct path match first (guarda de path traversal — igual a serve_fig_rmsd)
    direct = (figs_base / fname).resolve()
    if direct.exists() and str(direct).startswith(base_resolved):
        mime = "image/svg+xml" if direct.suffix == ".svg" else "image/png"
        with open(direct, "rb") as fh:
            return Response(fh.read(), mimetype=mime, headers={"Cache-Control":"max-age=300"})
    # Search in compound/session subdirs
    for compound_dir in sorted(figs_base.iterdir()) if figs_base.exists() else []:
        if not compound_dir.is_dir(): continue
        for sess_dir in sorted(compound_dir.iterdir(), reverse=True):
            if not sess_dir.is_dir(): continue
            candidate = (sess_dir / fname).resolve()
            if candidate.exists() and str(candidate).startswith(base_resolved):
                mime = "image/svg+xml" if candidate.suffix == ".svg" else "image/png"
                with open(candidate, "rb") as fh:
                    return Response(fh.read(), mimetype=mime, headers={"Cache-Control":"max-age=300"})
    return Response("not found", status=404)


@app.route('/fig_rmsd/<path:fname>')
def serve_fig_rmsd(fname):
    """Serve figuras de RMSD/HADDOCK (v2.1) de docking_project/figures/ relativo
    ao script (independente de NQ_DATA_DIR)."""
    base = SCRIPT_DIR / "docking_project" / "figures"
    target = (base / fname).resolve()
    # seguranca: impedir path traversal
    if not str(target).startswith(str(base.resolve())) or not target.exists():
        return Response("not found", status=404)
    mime = "image/svg+xml" if target.suffix == ".svg" else "image/png"
    with open(target, "rb") as fh:
        return Response(fh.read(), mimetype=mime, headers={"Cache-Control": "max-age=120"})


@app.route('/md_rmsd_tags')
def md_rmsd_tags_route():
    """Lista as tags de MD com analise pronta (md_<tag>_rmsd_*) — integra a aba
    RMSD & HADDOCK com os resultados gerados por md_simulation.py::analyze()."""
    mdir = SCRIPT_DIR / "docking_project" / "figures" / "md_clean"
    tags = []
    if mdir.exists():
        for p in sorted(mdir.glob("md_*_rmsd_summary.json")):
            t = p.name[len("md_"):-len("_rmsd_summary.json")]
            if t not in tags:
                tags.append(t)
        for p in sorted(mdir.glob("md_*_rmsd_rmsf.png")):   # figuras sem json
            t = p.name[len("md_"):-len("_rmsd_rmsf.png")]
            if t not in tags:
                tags.append(t)
    return jsonify({"tags": tags})


@app.route('/md_rmsd_summary/<tag>')
def md_rmsd_summary_route(tag):
    """Resumo numerico (JSON) da analise de MD para uma tag (proteina/ligante/RMSF)."""
    import json as _json
    p = (SCRIPT_DIR / "docking_project" / "figures" / "md_clean"
         / f"md_{tag}_rmsd_summary.json")
    if not p.exists():
        return jsonify({"found": False})
    try:
        d = _json.loads(p.read_text(encoding="utf-8"))
        d["found"] = True
        return jsonify(d)
    except Exception as e:
        return jsonify({"found": False, "error": str(e)})


@app.route('/binding_mode_status')
def binding_mode_status_route():
    """Consolida o estado MAIS ATUAL das analises de modo de ligacao da MD:
    SUBVERSIVO (FAD si-face) e INIBIDOR (fenda periferica), lendo os JSONs de
    binding-mode do kit fig13_3d_kit/redock. Informa se cada analise foi feita,
    a data (mtime), metricas-chave (LIG->FAD, LIG->ditiol, contatos Cys, RMSD),
    e a pre-organizacao covalente (subversivo). Usado pela aba RMSD & HADDOCK."""
    import json as _json, datetime as _dt
    # Diretorio do kit (override por env NQ_KIT_REDOCK; default no projeto do paper).
    kit = os.environ.get("NQ_KIT_REDOCK") or (
        r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-snc-nq2"
        r"\montagem_resultados\fig13_3d_kit\redock")
    kitp = Path(kit)

    def _newest(patterns):
        cands = []
        if kitp.exists():
            for pat in patterns:
                cands += list(kitp.glob(pat))
        return max(cands, key=lambda p: p.stat().st_mtime) if cands else None

    def _summarize(path):
        if not path or not path.exists():
            return {"done": False}
        try:
            d = _json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            return {"done": False, "error": str(e)}
        mt = _dt.datetime.fromtimestamp(path.stat().st_mtime)
        contacts = d.get("top_contacts", []) or []
        cys = [[n, o] for (n, o) in contacts if str(n).upper().startswith("CYS")]
        def _mean(k):
            v = d.get(k) or {}
            return v.get("mean")
        return {
            "done": True, "source": path.name,
            "date": mt.strftime("%Y-%m-%d %H:%M"),
            "n_frames": d.get("n_frames"),
            "lig_fad_mean": _mean("lig_fad_A"),
            "lig_dithiol_mean": _mean("lig_dithiol_A"),
            "lig_rmsd_mean": _mean("lig_rmsd_A"),
            "cys_contacts": cys,
            "top_contacts": contacts[:8],
        }

    sub = _summarize(_newest(["binding_mode_*subversive*.json"]))
    inh = _summarize(_newest(["binding_mode_*inhibitor*.json", "binding_mode_*inibidor*.json"]))

    # Pre-organizacao covalente (so faz sentido p/ o subversivo)
    cov = {}
    covp = kitp / "covalent_preorganization.json"
    if covp.exists():
        try:
            c = _json.loads(covp.read_text(encoding="utf-8"))
            cov = {"sg_c2_mean": (c.get("sg_c2_A") or {}).get("mean"),
                   "burgi_dunitz_mean": (c.get("burgi_dunitz_deg") or {}).get("mean"),
                   "cys_resi": c.get("cys52_resi"),
                   "date": _dt.datetime.fromtimestamp(covp.stat().st_mtime).strftime("%Y-%m-%d %H:%M")}
        except Exception:
            cov = {}

    return jsonify({
        "subversive": sub, "inhibitor": inh, "covalent": cov,
        "both_done": bool(sub.get("done") and inh.get("done")),
        "kit_dir": str(kitp),
    })


@app.route('/guide_pdf')
def guide_pdf_route():
    """Serve o Guia de Estudo (Docking + Dinamica Molecular) em PDF."""
    p = SCRIPT_DIR / "referencia" / "Guia_Docking_e_Dinamica_Molecular.pdf"
    if not p.exists():
        return Response("Guia PDF nao encontrado em referencia/", status=404)
    with open(p, "rb") as fh:
        return Response(fh.read(), mimetype="application/pdf",
                        headers={"Content-Disposition":
                                 "inline; filename=Guia_Docking_e_Dinamica_Molecular.pdf"})


@app.route('/run_md', methods=['POST'])
def run_md_route():
    # MD module aposentado em 2026-05-12; use /run_md_clean
    return jsonify({
        "log": "[deprecated] /run_md aposentado em 2026-05-12. Use /run_md_clean (botao 'Rodar MD Clean' na aba MD Simulation).",
        "results": {}, "text": "", "deprecated": True
    }), 410



# ═════════════════════════════════════════════════════════════════════════
#  MD CLEAN — endpoints para md_clean.py (subprocess CUDA, sem linha de comando)
# ═════════════════════════════════════════════════════════════════════════
_md_clean_proc = None  # subprocess.Popen handle (global, single-user dashboard)
_md_power_cap_to_revert = None  # W p/ reverter o power-cap da GPU quando o MD terminar


def _set_gpu_power_limit(watts):
    """Aplica nvidia-smi -pl <watts> (precisa admin). Retorna (ok, msg)."""
    import subprocess as _sp
    try:
        r = _sp.run(["nvidia-smi", "-pl", str(int(watts))],
                    capture_output=True, text=True, timeout=20)
        ok = (r.returncode == 0)
        msg = (r.stdout or "") + (r.stderr or "")
        return ok, msg.strip()
    except Exception as e:
        return False, str(e)

@app.route('/run_md_clean', methods=['POST'])
def run_md_clean_route():
    """Dispara md_clean.py como subprocess (nao-blocking) com env 'datascience'.
    Body JSON: {ns: float, platform: "CUDA"|"OpenCL"|"CPU"|"AUTO", tag: str}
    Retorna: {started: bool, pid: int, log_path: str}
    """
    global _md_clean_proc
    import subprocess, shutil
    data = request.json or {}
    ns       = float(data.get('ns', 5.0))
    platform = str(data.get('platform', 'CUDA')).strip()
    tag      = str(data.get('tag', 'TcTR')).strip().replace(' ', '_') or 'TcTR'

    # ja rodando?
    if _md_clean_proc is not None and _md_clean_proc.poll() is None:
        return jsonify({"started": False,
                        "error": f"MD Clean ja em execucao (PID {_md_clean_proc.pid})."}), 409

    script = SCRIPT_DIR / "md_simulation.py"
    if not script.exists():
        script = SCRIPT_DIR / "md_clean.py"  # fallback p/ shim de compatibilidade
    if not script.exists():
        return jsonify({"started": False, "error": f"md_simulation.py nao encontrado em {SCRIPT_DIR}"}), 500

    # Python do env datascience (mesmo padrao do venv_info.json)
    env_root = Path(r"C:\Users\luizd\miniforge3\envs\datascience")
    py_exe = env_root / "python.exe"
    if not py_exe.exists():
        # fallback: python atual (se dashboard ja roda dentro do env)
        py_exe = Path(sys.executable)
        env_root = Path(sys.prefix)

    target  = str(data.get('target', 'TcTR')).strip() or 'TcTR'
    control = str(data.get('control', 'none')).strip() or 'none'
    lang    = str(data.get('lang', 'pt')).strip()[:2] or 'pt'

    # FRESH como default (FAST exige complex_<TARGET>_initial.pdb pre-gerado,
    # que so existe para TcTR — e mesmo assim sem FAD: por isso o default mudou).
    fresh = bool(data.get('fresh', True))
    mode_flag = "--fresh" if fresh else "--fast"
    # Hipotese mecanistica: 'inhibitor' (pose periferica) ou 'subversive' (si-face do FAD)
    bind_mode = str(data.get('mode', 'inhibitor')).strip().lower()
    if bind_mode not in ('inhibitor', 'subversive'):
        bind_mode = 'inhibitor'
    cmd = [str(py_exe), str(script), mode_flag,
           "--ns", str(ns),
           "--tag", tag,
           "--platform", platform,
           "--target", target,
           "--mode", bind_mode]
    if platform == "AUTO":
        cmd.append("--allow-fallback")

    log_dir = SCRIPT_DIR / "docking_project" / "figures" / "md_clean"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "md_clean.log"

    # ── ENV: replicar `conda activate datascience` para o subprocess ─────────
    # OpenMM CUDA platform precisa ver cudart64_*.dll que ficam em
    # <env>\Library\bin. Se o dashboard rodar fora do env ativado, o subprocess
    # nao herdaria essas DLLs e cairia para CPU/OpenCL silenciosamente.
    sub_env = os.environ.copy()
    cuda_paths = [
        str(env_root),                              # python.exe
        str(env_root / "Library" / "mingw-w64" / "bin"),
        str(env_root / "Library" / "usr" / "bin"),
        str(env_root / "Library" / "bin"),          # ★ cudart64_*.dll, cufft, etc.
        str(env_root / "Scripts"),
        str(env_root / "bin"),
    ]
    existing_path = sub_env.get("PATH", "")
    sub_env["PATH"] = ";".join([p for p in cuda_paths if Path(p).exists()] + [existing_path])
    sub_env["CONDA_PREFIX"] = str(env_root)
    sub_env["CONDA_DEFAULT_ENV"] = "datascience"
    # Garante que OpenMM ache plugins CUDA (alguns wheels usam OPENMM_PLUGIN_DIR)
    plugin_dir = env_root / "Library" / "lib" / "plugins"
    if plugin_dir.exists():
        sub_env["OPENMM_PLUGIN_DIR"] = str(plugin_dir)

    # ── POWER-CAP da GPU (evita estourar o no-break / Kernel-Power 41) ───────
    # power_cap: "full" (padrao 170 W) ou "120"/"100" (reduzido). Precisa ADMIN
    # para o -pl funcionar; se falhar, segue mesmo assim (so loga o aviso).
    global _md_power_cap_to_revert
    power_cap = str(data.get('power_cap', 'full')).strip().lower()
    cap_note = ""
    if power_cap not in ('full', '', 'none'):
        try:
            watts = int(power_cap)
        except Exception:
            watts = 120
        ok, msg = _set_gpu_power_limit(watts)
        if ok:
            _md_power_cap_to_revert = 170   # reverte ao padrao quando terminar
            cap_note = f"[power-cap] GPU limitada a {watts} W (revertera p/ 170 W ao terminar)."
        else:
            cap_note = (f"[power-cap] FALHOU aplicar {watts} W (rode o dashboard como "
                        f"ADMINISTRADOR). Detalhe: {msg[:200]}")
    else:
        _md_power_cap_to_revert = None
        cap_note = "[power-cap] potencia TOTAL (sem limite extra)."
    print("[run_md_clean] " + cap_note)

    try:
        log_fh = open(log_path, "wb", buffering=0)
        if cap_note:
            log_fh.write((cap_note + "\n").encode("utf-8"))
        _md_clean_proc = subprocess.Popen(
            cmd, cwd=str(SCRIPT_DIR), env=sub_env,
            stdout=log_fh, stderr=subprocess.STDOUT,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        print(f"[run_md_clean] PID {_md_clean_proc.pid} cmd={' '.join(cmd)}")
        print(f"[run_md_clean] env CONDA_PREFIX={sub_env['CONDA_PREFIX']}")
        print(f"[run_md_clean] env OPENMM_PLUGIN_DIR={sub_env.get('OPENMM_PLUGIN_DIR','(none)')}")
    except Exception as e:
        return jsonify({"started": False, "error": f"Popen falhou: {e}"}), 500

    return jsonify({"started": True,
                    "pid": _md_clean_proc.pid,
                    "log_path": str(log_path),
                    "cmd": " ".join(cmd)})


@app.route('/md_clean_progress')
def md_clean_progress_route():
    """Le tail do md_clean.log. Indica termino via _md_clean_proc.poll()."""
    global _md_clean_proc
    n_lines = int(request.args.get('lines', 200))
    log_path = SCRIPT_DIR / "docking_project" / "figures" / "md_clean" / "md_clean.log"

    txt = ""
    if log_path.exists():
        try:
            with open(log_path, "rb") as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(0, size - 64 * 1024))   # ultimos ~64 KB
                raw = fh.read().decode("utf-8", errors="replace")
            lines = raw.splitlines()
            txt = "\n".join(lines[-n_lines:])
        except Exception as e:
            txt = f"[md_clean_progress] erro lendo log: {e}"

    done, exit_code, pid = True, None, None
    if _md_clean_proc is not None:
        pid = _md_clean_proc.pid
        rc = _md_clean_proc.poll()
        if rc is None:
            done = False
        else:
            exit_code = rc

    # Reverte o power-cap da GPU quando o MD termina (so uma vez).
    global _md_power_cap_to_revert
    if done and _md_power_cap_to_revert:
        ok, _ = _set_gpu_power_limit(_md_power_cap_to_revert)
        print(f"[md_clean_progress] power-cap revertido p/ {_md_power_cap_to_revert} W (ok={ok}).")
        _md_power_cap_to_revert = None

    return jsonify({"log": txt, "done": done, "exit_code": exit_code, "pid": pid})


@app.route('/md_clean_metrics')
def md_clean_metrics_route():
    """Serie temporal AO VIVO da MD, lida do CSV do StateDataReporter (md_<tag>.csv).
    Alimenta os graficos que atualizam durante a simulacao no dashboard. Funciona
    tanto para runs iniciados pelo dashboard quanto pela linha de comando."""
    import csv as _csv, re as _re
    tag = str(request.args.get('tag', 'TcTR_run01')).strip().replace(' ', '_')
    # seguranca: whitelist neutraliza '/', '\\' e '..' (path traversal via query)
    tag = _re.sub(r'[^A-Za-z0-9_.-]', '_', tag) or 'TcTR_run01'
    try:
        max_pts = int(request.args.get('max', 240))
    except Exception:
        max_pts = 240
    csv_path = (SCRIPT_DIR / "docking_project" / "figures" / "md_clean" / f"md_{tag}.csv")
    if not csv_path.exists():
        return jsonify({"found": False, "tag": tag})

    def _num(s):
        try:
            return float(s)
        except (TypeError, ValueError):
            return None

    t_ns, temp, pe, te = [], [], [], []
    last = {}
    try:
        with open(csv_path, "r", encoding="utf-8", errors="replace") as fh:
            rows = list(_csv.reader(fh))
        if not rows:
            return jsonify({"found": False, "tag": tag})
        header = [c.lstrip('#').strip().strip('"') for c in rows[0]]

        def colidx(*names):
            for n in names:
                for i, h in enumerate(header):
                    if h.lower().startswith(n.lower()):
                        return i
            return None

        i_time  = colidx("Time (ps)", "Time")
        i_temp  = colidx("Temperature")
        i_pe    = colidx("Potential")
        i_te    = colidx("Total Energy", "Total")
        i_prog  = colidx("Progress")
        i_speed = colidx("Speed")
        i_rem   = colidx("Time Remaining")

        def cell(row, i):
            return row[i] if (i is not None and i < len(row)) else None

        for row in rows[1:]:
            if not row:
                continue
            tv = _num(cell(row, i_time))
            if tv is None:
                continue
            t_ns.append(round(tv / 1000.0, 4))   # ps -> ns
            temp.append(_num(cell(row, i_temp)))
            pe.append(_num(cell(row, i_pe)))
            te.append(_num(cell(row, i_te)))
            last = {
                "progress": cell(row, i_prog),
                "speed": _num(cell(row, i_speed)),
                "time_remaining": cell(row, i_rem),
                "temp": temp[-1], "te": te[-1], "t_ns": t_ns[-1],
            }
    except Exception as e:
        return jsonify({"found": False, "tag": tag, "error": str(e)})

    n = len(t_ns)
    if n > max_pts and max_pts > 1:                # downsample preservando o ultimo ponto
        stepf = n / float(max_pts)
        keep = sorted(set([int(k * stepf) for k in range(max_pts)] + [n - 1]))
        t_ns = [t_ns[k] for k in keep]
        temp = [temp[k] for k in keep]
        pe   = [pe[k]   for k in keep]
        te   = [te[k]   for k in keep]
    running = bool(_md_clean_proc is not None and _md_clean_proc.poll() is None)
    return jsonify({"found": True, "tag": tag, "n": n, "running": running,
                    "t_ns": t_ns, "temp": temp, "pe": pe, "te": te, "last": last})


@app.route('/md_clean_stop', methods=['POST'])
def md_clean_stop_route():
    global _md_clean_proc
    if _md_clean_proc is None or _md_clean_proc.poll() is not None:
        return jsonify({"stopped": False, "message": "nenhum processo MD Clean ativo"})
    try:
        _md_clean_proc.terminate()
        try:
            _md_clean_proc.wait(timeout=5)
        except Exception:
            _md_clean_proc.kill()
        return jsonify({"stopped": True, "message": f"PID {_md_clean_proc.pid} encerrado"})
    except Exception as e:
        return jsonify({"stopped": False, "message": f"erro: {e}"}), 500


@app.route('/md_clean_script', methods=['POST'])
def md_clean_script_route():
    """Devolve o md_clean.py com cabecalho personalizado para o usuario copiar/baixar."""
    data = request.json or {}
    ns       = float(data.get('ns', 5.0))
    platform = str(data.get('platform','CUDA')).strip()
    tag      = str(data.get('tag','TcTR')).strip().replace(' ','_') or 'TcTR'
    target   = str(data.get('target','TcTR')).strip() or 'TcTR'
    lang     = str(data.get('lang','pt')).strip()[:2] or 'pt'

    src = SCRIPT_DIR / "md_simulation.py"
    if not src.exists():
        src = SCRIPT_DIR / "md_clean.py"  # fallback p/ shim
    body = src.read_text(encoding='utf-8') if src.exists() else "# md_simulation.py nao encontrado\n"

    if lang == 'pt':
        header = (
            "# ========================================================\n"
            f"# Script OpenMM gerado pelo dashboard NQ-Docking\n"
            f"# Alvo:        {target}\n"
            f"# Tempo prod:  {ns} ns\n"
            f"# Plataforma:  {platform}\n"
            f"# Tag:         {tag}\n"
            f"# Idioma:      {lang}\n"
            f"# ========================================================\n"
            f"# Para reproduzir esta simulacao:\n"
            f"#   conda activate datascience\n"
            f"#   python md_simulation.py --fast --ns {ns} --tag {tag} --platform {platform} --target {target}\n"
            f"# ========================================================\n\n"
        )
    else:
        header = (
            "# ========================================================\n"
            f"# OpenMM script generated by NQ-Docking dashboard\n"
            f"# Target:      {target}\n"
            f"# Production:  {ns} ns\n"
            f"# Platform:    {platform}\n"
            f"# Tag:         {tag}\n"
            f"# Language:    {lang}\n"
            f"# ========================================================\n"
            f"# To reproduce:\n"
            f"#   conda activate datascience\n"
            f"#   python md_simulation.py --fast --ns {ns} --tag {tag} --platform {platform} --target {target}\n"
            f"# ========================================================\n\n"
        )
    from flask import Response
    return Response(header + body, mimetype="text/plain")


@app.route('/md_progress')
def md_progress_route():
    """Retorna o conteúdo atual do arquivo de progresso live do MD.
    Usado pelo dashboard para polling durante simulação em tempo real.
    Parâmetros: ?target=TcTR&session=...
    """
    target      = request.args.get('target', 'TcTR')
    session_tag = request.args.get('session', None)
    compound    = request.args.get('compound', '2-SCN-NQ')

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        # Tentar pegar a sessão mais recente
        results_dir = BASE_DIR / "results"
        if results_dir.exists():
            subdirs = sorted(results_dir.iterdir(), reverse=True)
            session_tag = subdirs[0].name if subdirs else make_session_tag()
        else:
            session_tag = make_session_tag()

    dirs = get_session_dirs(BASE_DIR, session_tag, compound)
    prog_file = dirs["figures"] / f"md_progress_{target}.log"

    if not prog_file.exists():
        return jsonify({"lines": [], "running": False})

    try:
        content = prog_file.read_text(encoding="utf-8", errors="replace")
        lines   = [l for l in content.splitlines() if l.strip()]
        # Detecta se ainda está rodando: lock file presente
        lock_file = dirs["figures"] / f"md_{target}.lock"
        running   = lock_file.exists()
        return jsonify({"lines": lines, "running": running, "target": target})
    except Exception as e:
        return jsonify({"lines": [f"[ERRO] {e}"], "running": False})


@app.route('/md_status')
def md_status_route():
    """Return OpenMM/CUDA diagnostic info for the dashboard."""
    return jsonify({
        "deprecated": True,
        "message": "md_module aposentado — use /md_clean_progress",
        "python_exe": sys.executable,
    })


@app.route('/check_docking_for_md', methods=['POST'])
def check_docking_for_md_route():
    """Verifica se cada alvo selecionado tem result_<TID>.pdbqt (pose docked do Vina).

    Retorna:
      {
        "missing":  [{"target": "TcTR", "reason": "no_docked_pose"}, ...],
        "ok":       [{"target": "LiTR", "path": "..."}],
        "compound": "2-SCN-NQ"
      }

    Sem pose docked, MD em modo 'real' produz resultados FISICAMENTE INVALIDOS
    (ligante posicionado arbitrariamente fora do bolsao). Este gate bloqueia
    a execucao da MD ate que o docking seja realizado.
    """
    data = request.json or {}
    targets = data.get('targets', []) or []
    compound = _clean_compound(data.get('compound'))

    _ensure_script_path()
    from compound_db import _safe_name
    safe_c = _safe_name(compound)
    script_dir = SCRIPT_DIR

    # Mesmo conjunto de raizes usado em md_module.run_md_analysis
    base = Path("docking_project")
    base_abs = base.resolve() if base.is_absolute() else (script_dir / base).resolve()
    roots = [
        base_abs,
        script_dir / "docking_project",
        script_dir,
        Path(r"E:\data\docking"),
        Path(r"E:\data\docking") / "docking_project",
    ]
    seen, roots_unique = set(), []
    for r in roots:
        rk = str(r).lower()
        if rk not in seen:
            seen.add(rk); roots_unique.append(r)

    missing = []
    ok_list = []

    def find_docked(tid):
        for root in roots_unique:
            if not root.exists():
                continue
            # 1) docking/<compound>/<session>/result_<TID>.pdbqt
            dock_base = root / "docking" / safe_c
            if dock_base.exists():
                try:
                    sessions = sorted(
                        [s for s in dock_base.iterdir() if s.is_dir()], reverse=True
                    )
                except Exception:
                    sessions = []
                for sess in sessions:
                    cand = sess / f"result_{tid}.pdbqt"
                    if cand.exists() and cand.stat().st_size > 100:
                        return ("docked_session", str(cand))
            # 2) docking/result_<TID>.pdbqt (flat)
            flat = root / "docking" / f"result_{tid}.pdbqt"
            if flat.exists() and flat.stat().st_size > 100:
                return ("docked_flat", str(flat))
        return (None, None)

    for tid in targets:
        source, path = find_docked(tid)
        if source is None:
            missing.append({"target": tid, "reason": "no_docked_pose"})
        else:
            ok_list.append({"target": tid, "source": source, "path": path})

    return jsonify({
        "missing": missing,
        "ok": ok_list,
        "compound": compound,
        "all_ok": len(missing) == 0,
    })


@app.route('/run_docking', methods=['POST'])
def run_docking_route():
    """Roda docking_v2.py para os alvos especificados.

    Payload JSON:
      {
        "targets": ["TcTR", "LiTR", ...],        # ou "all" / "parasite"
        "exhaustiveness": 32,                      # padrao 32
        "compound": "2-SCN-NQ"                     # informativo
      }

    Retorna:
      {
        "log":        "...stdout/stderr completo...",
        "returncode": 0 | non-zero,
        "results":    {tid: {best_dG, ki_uM, ...}}  (se docking_results.json existe)
      }
    """
    data = request.json or {}
    targets = data.get('targets', None)
    exhaustiveness = int(data.get('exhaustiveness', 32))
    compound = _clean_compound(data.get('compound'))

    # Construir comando
    py_exe = sys.executable
    script = SCRIPT_DIR / "docking_v2.py"
    cmd = [py_exe, str(script), "--exhaustiveness", str(exhaustiveness)]

    if isinstance(targets, str):
        if targets.lower() == "all":
            cmd.append("--all")
        elif targets.lower() == "parasite":
            cmd.append("--parasite")
        else:
            cmd += ["--targets", targets]
    elif isinstance(targets, list) and targets:
        cmd += ["--targets", ",".join(targets)]
    else:
        return jsonify({"error": "targets nao especificado"}), 400

    log_buf = io.StringIO()
    log_buf.write(f"[INFO] Iniciando docking de {compound}\n")
    log_buf.write(f"[INFO] Comando: {' '.join(cmd)}\n")
    log_buf.write(f"[INFO] Tempo estimado: ~5-10 min por alvo (exhaustiveness={exhaustiveness})\n\n")

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(SCRIPT_DIR),
            timeout=10800,  # 3 horas max
        )
        log_buf.write(proc.stdout or "")
        if proc.stderr:
            log_buf.write("\n--- STDERR ---\n")
            log_buf.write(proc.stderr)
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        log_buf.write("\n[ERRO] Timeout (>3h) — docking abortado.\n")
        rc = -1
    except Exception as e:
        log_buf.write(f"\n[ERRO] Falha ao executar: {e}\n")
        rc = -2

    # Tentar ler resultados JSON se existirem.
    # Antes buscava um caminho fixo relativo ao SCRIPT, ignorando o proprio
    # NQ_DATA_DIR do dashboard: quem apontasse os dados para outro disco rodava
    # o docking e via "sem resultado". results_io conhece as convencoes todas.
    results = {}
    try:
        import results_io as _rio
        results = _rio.read_results(_default_campaign(), compound) or {}
    except Exception as e:
        app.logger.warning("docking: falha ao ler o resultado agregado: %s", e)
        try:
            results_json = SCRIPT_DIR / "docking_project" / "results" / "docking_results.json"
            if results_json.exists():
                results = json.loads(results_json.read_text(encoding="utf-8"))
        except Exception as e2:
            app.logger.warning("docking: fallback legado tambem falhou: %s", e2)

    return jsonify({
        "log": log_buf.getvalue(),
        "returncode": rc,
        "results": results,
        "compound": compound,
    })


@app.route('/help_doc/<path:key>')
def help_doc_route(key):
    """Serve markdown doc content for the Ajuda tab (docs/ consolidado)."""
    import pathlib
    doc_map = {
        'readme':     'README.md',
        'protocolo':  'PROTOCOLO_CIENTIFICO.md',
        'drugdesign': 'DRUG_DESIGN_MANUAL.md',
        'proxy':      'PROXY_REVERSO_CLOUDFLARE.md',
    }
    fname = doc_map.get(key, key)
    # seguranca: neutraliza path traversal ('..', separadores) — so o nome base
    fname = os.path.basename(fname.replace('\\', '/'))
    if not fname or not fname.lower().endswith('.md'):
        fname = (fname + '.md') if fname else 'README.md'
    docs_dir = SCRIPT_DIR / 'docs'
    candidates = [docs_dir / fname, SCRIPT_DIR / fname,
                  SCRIPT_DIR.parent / fname, pathlib.Path('docs') / fname,
                  pathlib.Path('.') / fname]
    for path in candidates:
        try:
            if path.exists():
                return jsonify({"raw": path.read_text(encoding='utf-8'), "found": True})
        except Exception:
            continue
    return jsonify({"error": f"Arquivo não encontrado: {fname}. "
                              f"Verifique a pasta docs/ ao lado do dashboard_v2.py.", "raw": ""})


@app.route('/help_list')
def help_list_route():
    """Lista todos os .md consolidados em docs/ (+ README.md da raiz) para a aba Ajuda."""
    docs_dir = SCRIPT_DIR / 'docs'
    items = []
    readme = SCRIPT_DIR / 'README.md'
    if readme.exists():
        items.append({"file": "README.md", "title": "README"})
    if docs_dir.exists():
        for p in sorted(docs_dir.glob('*.md')):
            title = p.stem
            try:
                for line in p.read_text(encoding='utf-8', errors='replace').splitlines():
                    s = line.strip()
                    if s.startswith('#'):
                        title = s.lstrip('#').strip() or p.stem
                        break
            except Exception as e:
                app.logger.debug("help docs: falha ao ler titulo de %s: %s", p.name, e)
            items.append({"file": p.name, "title": title})
    return jsonify({"docs": items})


# ── Engenharia DUAL (subversivo + inibidor) — pipeline local, sem tokens ─────────
_dual_design_proc = None


@app.route('/run_dual_design', methods=['POST'])
def run_dual_design_route():
    """Dispara dual_design_pipeline.py em background (RDKit local; ollama opcional)."""
    global _dual_design_proc
    import subprocess
    if _dual_design_proc is not None and _dual_design_proc.poll() is None:
        return jsonify({"started": False, "error": f"Ja rodando (PID {_dual_design_proc.pid})."}), 409
    script = SCRIPT_DIR / "dual_design_pipeline.py"
    if not script.exists():
        return jsonify({"started": False, "error": "dual_design_pipeline.py nao encontrado."}), 500
    py = Path(r"C:\Users\luizd\miniforge3\envs\datascience\python.exe")
    if not py.exists():
        py = Path(sys.executable)
    cmd = [str(py), str(script)]
    if (request.json or {}).get('ollama'):
        cmd.append("--ollama")
    try:
        _dual_design_proc = subprocess.Popen(cmd, cwd=str(SCRIPT_DIR),
                                             creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    except Exception as e:
        return jsonify({"started": False, "error": f"Popen falhou: {e}"}), 500
    return jsonify({"started": True, "pid": _dual_design_proc.pid})


@app.route('/dual_design_status')
def dual_design_status_route():
    global _dual_design_proc
    running = _dual_design_proc is not None and _dual_design_proc.poll() is None
    rep = SCRIPT_DIR / "docking_project" / "figures" / "dual_design" / "dual_design_report.html"
    return jsonify({"running": running, "report_ready": rep.exists()})


@app.route('/dual_design_report')
def dual_design_report_route():
    rep = SCRIPT_DIR / "docking_project" / "figures" / "dual_design" / "dual_design_report.html"
    if not rep.exists():
        return ("<p style='font-family:sans-serif;padding:20px'>Relatorio ainda nao gerado. "
                "Clique em <b>Gerar 5 estruturas duais</b>.</p>")
    return rep.read_text(encoding="utf-8")


@app.route('/run_drug_design', methods=['POST'])
def run_drug_design_route():
    data           = request.get_json(silent=True) or {}
    approach       = data.get('approach', 'fbdd')
    parent_smiles  = data.get('parent_smiles', 'O=C1C=C(SC#N)C(=O)c2ccccc21')
    parent_name    = data.get('parent_name', '2-SCN-NQ')
    lang           = data.get('lang', 'en')
    compound_name  = _clean_compound(data.get('compound_name'))
    session_tag    = data.get('session_tag') or None

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        session_tag = make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)

    log_buf = io.StringIO(); old_out = sys.stdout; sys.stdout = log_buf
    results = {}
    try:
        from importlib import reload
        import drug_design_module; reload(drug_design_module)

        kwargs = dict(
            parent_smiles = parent_smiles,
            parent_name   = parent_name,
            approach      = approach,
            output_dir    = str(dirs["figures"]),
            lang          = lang,
        )
        if approach == 'fbdd':
            kwargs['rgroup_categories'] = data.get('categories') or None
            kwargs['max_variants']      = int(data.get('max_variants', 20))
        elif approach == 'hybrid':
            kwargs['partner_smiles'] = data.get('partner_smiles', '')
            kwargs['partner_name']   = data.get('partner_name', 'Partner')
            kwargs['max_variants']   = int(data.get('max_variants', 10))
        elif approach == 'ga':
            kwargs['seed_smiles']  = data.get('seed_smiles', [parent_smiles])
            kwargs['seed_names']   = data.get('seed_names',  [parent_name])
            kwargs['generations']  = int(data.get('generations', 15))
            kwargs['pop_size']     = int(data.get('pop_size', 20))

        results = drug_design_module.run_drug_design(**kwargs)

        # Save text
        if results.get('text'):
            (dirs["texts"] / f"drug_design_{approach}_{lang}.txt").write_text(
                results['text'], encoding='utf-8')
    except Exception as e:
        import traceback
        log_buf.write(f"\n[ERRO DrugDesign] {e}\n{traceback.format_exc()}")
    finally:
        sys.stdout = old_out

    # Serialise — strip non-JSON-safe objects
    def _safe(obj):
        if isinstance(obj, dict):
            return {k: _safe(v) for k, v in obj.items()
                    if k not in ('admet',)}
        if isinstance(obj, list):
            return [_safe(i) for i in obj]
        try:
            json.dumps(obj)
            return obj
        except Exception:
            return str(obj)

    return jsonify({"log": log_buf.getvalue(), "session_tag": session_tag,
                    **_safe(results)})


@app.route('/run_cyclodextrin', methods=['POST'])
def run_cyclodextrin_route():
    data            = request.get_json(silent=True) or {}
    lang            = data.get('lang','en')
    compound_name   = _clean_compound(data.get('compound_name'))
    compound_smiles = data.get('compound_smiles', None)
    session_tag     = data.get('session_tag') or None

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag, load_compound
    if session_tag is None:
        session_tag = make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)

    log_buf = io.StringIO(); old_out = sys.stdout; sys.stdout = log_buf
    cd_data = {}; text_out = ""
    try:
        from importlib import reload
        import cyclodex_module; reload(cyclodex_module)
        cd_data = cyclodex_module.run_cyclodextrin(
            str(dirs["figures"]), lang,
            compound_name=compound_name,
            compound_abbr=compound_name,
            ligand_smiles=compound_smiles,
        )
        text_out = cyclodex_module.generate_cd_text(cd_data, lang=lang)
        (dirs["texts"] / f"cyclodextrin_{lang}.txt").write_text(text_out, encoding='utf-8')
    except Exception as e:
        import traceback; log_buf.write(f"\n[ERRO CD] {e}\n{traceback.format_exc()}")
    finally:
        sys.stdout = old_out
    return jsonify({"log": log_buf.getvalue(), "results": cd_data, "text": text_out,
                    "session_tag": session_tag})


@app.route('/run_toxicity', methods=['POST'])
def run_toxicity_route():
    data            = request.get_json(silent=True) or {}
    smiles          = data.get('smiles', "O=C1C=C(SC#N)C(=O)c2ccccc21")
    lang            = data.get('lang', 'en')
    compound_name   = _clean_compound(data.get('compound_name'))
    session_tag     = data.get('session_tag') or None
    exp_data        = data.get('exp_data') or {}   # NEW: experimental values from UI

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        session_tag = make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)

    log_buf = io.StringIO(); old_out = sys.stdout; sys.stdout = log_buf
    tox_data = {}; text_out = ""
    try:
        from importlib import reload
        import toxicity_module; reload(toxicity_module)
        tox_data = toxicity_module.run_toxicity(
            smiles, str(dirs["figures"]), lang,
            compound_name=compound_name,
            exp_data=exp_data if exp_data else None,   # pass experimental values
        )
        text_out = toxicity_module.generate_toxicity_text(tox_data, lang=lang)
        (dirs["texts"] / f"toxicity_{lang}.txt").write_text(text_out, encoding='utf-8')
        # Save experimental data alongside results
        if exp_data:
            exp_path = dirs["results"] / "exp_data.json"
            exp_path.parent.mkdir(parents=True, exist_ok=True)
            with open(exp_path, 'w') as fh:
                json.dump({"compound": compound_name, "session": session_tag,
                           "experimental": exp_data}, fh, indent=2)
    except Exception as e:
        import traceback; log_buf.write(f"\n[ERRO TOX] {e}\n{traceback.format_exc()}")
    finally:
        sys.stdout = old_out
    return jsonify({"log": log_buf.getvalue(), "results": tox_data, "text": text_out,
                    "session_tag": session_tag})


@app.route('/run_all_figures', methods=['POST'])
def run_all_figures():
    data            = request.get_json(silent=True) or {}
    lang            = data.get('lang','en')
    smiles          = data.get('smiles',"O=C1C=C(SC#N)C(=O)c2ccccc21")
    compound_name   = _clean_compound(data.get('compound_name'))
    session_tag     = data.get('session_tag') or None

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        session_tag = make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound_name)

    log_buf = io.StringIO(); old_out = sys.stdout; sys.stdout = log_buf
    generated = []
    for mod_name, func_name, kwargs in [
        ("admet_module",    "run_admet",        {"smiles":smiles,"output_dir":str(dirs["figures"]),"lang":lang}),
        ("md_module",       "run_md_analysis",  {"output_dir":str(dirs["figures"]),"lang":lang}),
        ("cyclodex_module", "run_cyclodextrin", {"output_dir":str(dirs["figures"]),"lang":lang}),
        ("toxicity_module", "run_toxicity",     {"smiles":smiles,"output_dir":str(dirs["figures"]),"lang":lang}),
    ]:
        try:
            from importlib import import_module, reload
            if mod_name == 'md_module':
                # md_module aposentado em 2026-05-12 — use MD Clean no dashboard.
                generated.append(f"[SKIP] {mod_name} (aposentado — use MD Clean)")
                continue
            mod = reload(import_module(mod_name))
            getattr(mod, func_name)(**kwargs)
            generated.append(f"[OK] {mod_name}")
        except Exception as e:
            generated.append(f"[ERRO] {mod_name}: {str(e)[:100]}")
    rf = dirs["results"] / "results.json"
    if rf.exists():
        try:
            from analysis import run_analysis
            run_analysis(str(rf), str(dirs["figures"]), lang)
            generated.append("[OK] analysis (docking figures)")
        except Exception as e:
            generated.append(f"[ERRO] analysis: {str(e)[:80]}")
    sys.stdout = old_out
    figs = [f.name for f in sorted(dirs["figures"].iterdir()) if f.suffix=='.png'] \
           if dirs["figures"].exists() else []
    return jsonify({"generated": generated, "figures": figs, "session_tag": session_tag})


# ─── COMPOUND DATABASE ROUTES ─────────────────────────────────────

def _compounds_from_registry(campaign_id):
    """Adapta registry.list_ligands(campaign) -> formato que /compounds espera.

    Usado para campanhas fora do layout legado do compound_db (SUB01). Enriquece
    _has_pdbqt/_has_mol2/_dir a partir do disco quando os arquivos preparados de
    fato existem (ex.: latrunculina_A do SUB04); caso contrario marca como pendente
    -- honesto (fail-visible), nunca inventa arquivo inexistente. Mapeia o schema
    2.0 do registry (mw_g_mol, topology.expected_active_torsions, prep.optimization)
    para os campos flat que o frontend renderiza (mw, rotatable_bonds, optimization).
    """
    import registry as _reg

    ldir_root = None
    try:
        import paths as _paths
        ldir_root = _paths.data_root(campaign_id, warn_legacy=False) / "ligand"
    except Exception:
        ldir_root = None

    out = []
    for lig in _reg.list_ligands(campaign_id):
        name = lig.get("name") or lig.get("display_name") or lig.get("id")
        topo = lig.get("topology") or {}
        prep = lig.get("prep") or {}
        prov = lig.get("provenance") or {}
        has_pdbqt = has_mol2 = False
        found_dir = None
        if ldir_root is not None:
            for cand in (name, lig.get("id"), lig.get("name")):
                if not cand:
                    continue
                try:
                    d = ldir_root / str(cand)
                    if d.exists():
                        found_dir = str(d)
                        has_pdbqt = any(d.glob("*.pdbqt"))
                        has_mol2 = any(d.glob("*.mol2"))
                        break
                except Exception:
                    pass
        out.append({
            "name": name,
            "display_name": lig.get("display_name"),
            "smiles": lig.get("smiles"),
            "formula": lig.get("formula"),
            "mw": lig.get("mw_g_mol"),
            "heavy_atoms": lig.get("heavy_atoms"),
            "rotatable_bonds": topo.get("expected_active_torsions"),
            "charge": lig.get("charge"),
            "optimization": prep.get("optimization"),
            "date_added": prov.get("date_added"),
            "notes": lig.get("notes"),
            "class": lig.get("class"),
            "_dir": found_dir,
            "_has_pdbqt": has_pdbqt,
            "_has_mol2": has_mol2,
            "_source": "registry",
        })
    return out


@app.route('/compounds')
def list_compounds_route():
    _ensure_script_path()
    campaign = request.args.get('campaign') or _default_campaign()
    try:
        if campaign == _default_campaign():
            # Layout legado (SUB01): fonte original preservada integralmente,
            # inclusive _has_pdbqt/_has_mol2/_dir reais e ensure_default_compound().
            from compound_db import list_compounds, ensure_default_compound
            ensure_default_compound()
            compounds = list_compounds()
        else:
            # Demais campanhas (ex.: SUB04 -- drogas do citoesqueleto): ligantes
            # declarados no registry, que ja e a fonte de verdade multi-campanha.
            compounds = _compounds_from_registry(campaign)
        return jsonify({"compounds": compounds, "campaign": campaign})
    except Exception as e:
        return jsonify({"compounds": [], "campaign": campaign, "error": str(e)})


@app.route('/populate_library', methods=['POST'])
def populate_library_route():
    """Popula o banco com a biblioteca built-in de NQs antiparasitárias."""
    _ensure_script_path()
    try:
        from compound_db import populate_library, clean_orphan_dirs
        removed = clean_orphan_dirs()
        results = populate_library(overwrite=False)
        return jsonify({"results": results, "removed_orphans": removed, "n": len(results)})
    except Exception as e:
        import traceback
        return jsonify({"error": str(e), "traceback": traceback.format_exc()[:600]})


@app.route('/add_compound', methods=['POST'])
def add_compound_route():
    data = request.get_json(silent=True) or {}
    name         = str(data.get('name','')).strip()
    smiles       = str(data.get('smiles','')).strip()
    formula      = data.get('formula','')
    display_name = data.get('display_name','')
    notes        = data.get('notes','')
    # Bug-guard: JS às vezes passa objeto serializado como nome
    if not name or name.startswith('{') or name == '[object Object]':
        return jsonify({"error": "name inválido — informe uma string (ex: '2-Cl-NQ')"}), 400
    if not smiles:
        return jsonify({"error": "smiles é obrigatório"}), 400
    _ensure_script_path()
    try:
        from compound_db import add_compound
        meta = add_compound(name, smiles, formula=formula,
                            display_name=display_name, notes=notes,
                            prepare_3d=True)
        return jsonify(meta)
    except Exception as e:
        import traceback
        return jsonify({"error": str(e), "traceback": traceback.format_exc()[:500]})


@app.route('/validate_smiles', methods=['POST'])
def validate_smiles():
    data   = request.get_json(silent=True) or {}
    smiles = data.get('smiles','')
    _ensure_script_path()
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"valid": False, "error": "RDKit could not parse SMILES"})
        return jsonify({
            "valid":       True,
            "formula":     rdMolDescriptors.CalcMolFormula(mol),
            "mw":          round(Descriptors.MolWt(mol), 3),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "hbd":         rdMolDescriptors.CalcNumHBD(mol),
            "hba":         rdMolDescriptors.CalcNumHBA(mol),
        })
    except ImportError:
        # RDKit not available — basic validation
        if len(smiles) < 3:
            return jsonify({"valid": False, "error": "SMILES too short"})
        return jsonify({"valid": True, "formula": "", "mw": 0, "heavy_atoms": 0,
                        "note": "RDKit not installed — full validation unavailable"})
    except Exception as e:
        return jsonify({"valid": False, "error": str(e)})


# ─── SESSION ROUTES ────────────────────────────────────────────────

@app.route('/sessions')
def list_sessions_route():
    compound = request.args.get('compound', '2-SCN-NQ')
    _ensure_script_path()
    try:
        from compound_db import list_sessions
        sessions = list_sessions(BASE_DIR, compound)
        # Annotate with actual directory existence
        for s in sessions:
            s["figures_exists"] = Path(s["figures_dir"]).exists() if s.get("figures_dir") else False
            s["docking_exists"] = s["figures_exists"]  # alias — same dir
            s["texts_exists"]   = Path(s["texts_dir"]).exists()   if s.get("texts_dir")   else False
        return jsonify({"sessions": sessions, "compound": compound})
    except Exception as e:
        return jsonify({"sessions": [], "error": str(e)})


@app.route('/session_results')
def session_results_route():
    compound    = request.args.get('compound', '2-SCN-NQ')
    session_tag = request.args.get('session', '')
    _ensure_script_path()
    from compound_db import get_session_dirs
    dirs = get_session_dirs(BASE_DIR, session_tag, compound)
    rf = dirs["results"] / "results.json"
    if rf.exists():
        with open(rf) as f:
            d = json.load(f)
        return jsonify({"results": d.get("results", {}), "meta": d})
    return jsonify({"results": {}, "error": "No results.json for this session"})


@app.route('/generate_paper_text', methods=['POST'])
def generate_paper_text():
    data = request.get_json(silent=True) or {}
    lang = data.get('lang','en'); section = data.get('section','results_docking')
    style = data.get('style','acs'); tone = data.get('tone','scientific')
    results_data = {}; admet_data = {}
    # Derive results.json from the session/compound, not a fixed path.
    # dirs["results"] is per-compound (base/results/<safe_compound>); session_tag
    # is accepted for API symmetry but does not affect the results location.
    compound    = _clean_compound(data.get('compound'))
    session_tag = data.get('session', '')
    try:
        _ensure_script_path()
        from compound_db import get_session_dirs
        dirs = get_session_dirs(BASE_DIR, session_tag, compound)
        rf = dirs["results"] / "results.json"
        if not rf.exists():
            # Last-resort: pick the most recent results.json across all compounds.
            cand = sorted((BASE_DIR/"results").glob("*/results.json"),
                          key=lambda p: p.stat().st_mtime, reverse=True) \
                   if (BASE_DIR/"results").exists() else []
            rf = cand[0] if cand else rf
        if rf.exists():
            with open(rf) as f: results_data = json.load(f).get("results",{})
        else:
            app.logger.warning("generate_paper_text: no results.json for compound=%r session=%r (looked in %s); "
                               "using hardcoded fallback affinities", compound, session_tag, rf)
    except Exception as e:
        app.logger.warning("generate_paper_text: failed to load session results.json (%s); using fallback", e)
    try:
        af = BASE_DIR/"results"/"admet_profile.json"
        if af.exists():
            with open(af) as f: admet_data = json.load(f)
    except Exception as e:
        app.logger.warning("generate_paper_text: falha ao ler %s: %s", af, e)
    if not admet_data:
        admet_data = {"MW":215.23,"LogP":2.41,"HBD":0,"HBA":3,"TPSA":67.72,"RotBonds":2,
                      "QED":0.61,"GI_absorption":"High","Caco2_Papp_nm_s":28.5,
                      "LD50_mouse_mg_kg":185,"LD50_human_pred":310,"Toxicity_class":3,
                      "CD_Ka_M1":1250,"Lipinski_pass":True}
    section_map = {
        "results_docking":"Results - Molecular Docking",
        "results_admet":"Results - In Silico ADMET Profile",
        "results_md":"Results - Molecular Dynamics Simulation",
        "results_cd":"Results - Beta-Cyclodextrin Complexation",
        "results_tox":"Results - Human Toxicity Prediction",
        "methods_docking":"Methods - Molecular Docking",
        "methods_admet":"Methods - ADMET Analysis",
        "discussion":"Discussion",
        "abstract":"Abstract / Resumo",
    }
    tone_map = {"scientific":"formal scientific prose","concise":"concise under 250 words","comprehensive":"comprehensive and detailed"}
    docking_summary = "TcTR: -9.2 kcal/mol, LiTR: -8.8 kcal/mol, hGR: -6.1 kcal/mol (SI~15)"
    if results_data:
        top3 = sorted(results_data.items(), key=lambda x: x[1].get('best_dG',0))[:3]
        docking_summary = ", ".join([f"{t}: dG={d.get('best_dG',0):.1f} kcal/mol" for t,d in top3])
    lang_name = "Portuguese (Brazilian)" if lang == "pt" else "English"
    prompt = (
        f"You are a scientific writer in medicinal chemistry and antiparasitic drug discovery.\n\n"
        f"Write the {section_map.get(section,section)} section for a manuscript about "
        f"2-Thiocyanato-1,4-naphthoquinone (2-SCN-NQ).\n\n"
        f"LANGUAGE: {lang_name}\nJOURNAL STYLE: {style}\nTONE: {tone_map.get(tone,'scientific')}\n\n"
        f"DOCKING: {docking_summary}\n"
        f"ADMET: MW={admet_data.get('MW',215.23)}, LogP={admet_data.get('LogP',2.41)}, "
        f"TPSA={admet_data.get('TPSA',67.72)}, QED={admet_data.get('QED',0.61)}, "
        f"Lipinski={'PASS' if admet_data.get('Lipinski_pass',True) else 'FAIL'}, "
        f"Caco-2={admet_data.get('Caco2_Papp_nm_s',28.5)} nm/s, "
        f"mouse LD50={admet_data.get('LD50_mouse_mg_kg',185)} mg/kg GHS Class {admet_data.get('Toxicity_class',3)}\n\n"
        f"KEY REFS: Trott & Olson 2010 (Vina); Lipinski 2001 (RO5); Bickerton 2012 (QED); "
        f"Banerjee 2018 (ProTox-II); Loftsson 2012 (beta-CD)\n\n"
        f"Write ONLY the section text. No markdown, no headers. (Author, year) citation style. "
        f"Length: {'150-250 words' if tone=='concise' else '300-500 words'}."
    )
    api_key = os.environ.get('ANTHROPIC_API_KEY')
    if not api_key:
        return jsonify({"error": "ANTHROPIC_API_KEY nao configurada no ambiente. Defina a variavel de ambiente ANTHROPIC_API_KEY (sk-ant-...) e reinicie o servidor, ou use 'Texto pre-definido'."})
    try:
        import urllib.request
        req_data = json.dumps({
            "model": "claude-sonnet-4-20250514",
            "max_tokens": 1000,
            "messages": [{"role":"user","content":prompt}]
        }).encode()
        req = urllib.request.Request(
            "https://api.anthropic.com/v1/messages",
            data=req_data,
            headers={
                "x-api-key": api_key,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json"
            },
            method="POST"
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            resp_data = json.loads(resp.read())
            text = resp_data.get("content",[{}])[0].get("text","")
        texts_dir = BASE_DIR/"manuscript_texts"; texts_dir.mkdir(parents=True, exist_ok=True)
        (texts_dir/f"{section}_{lang}_ai.txt").write_text(text, encoding='utf-8')
        return jsonify({"text": text})
    except Exception as e:
        return jsonify({"error": f"API error: {str(e)[:200]}. Use 'Texto pre-definido'."})


@app.route('/static_text', methods=['POST'])
def static_text():
    data = request.get_json(silent=True) or {}
    lang = data.get('lang','en'); section = data.get('section','results_admet')
    texts_dir = BASE_DIR / "manuscript_texts"
    key_map = {
        "results_admet":   f"ADMET_results_{lang}.txt",
        "results_docking": f"results_text_{lang}.txt",
        "results_md":      f"MD_results_{lang}.txt",
        "results_cd":      f"CD_results_{lang}.txt",
        "results_tox":     f"Toxicity_results_{lang}.txt",
    }
    fname = key_map.get(section)
    if fname and (texts_dir/fname).exists():
        return jsonify({"text": (texts_dir/fname).read_text(encoding='utf-8')})
    _ensure_script_path()
    try:
        if "admet" in section:
            import admet_module
            return jsonify({"text": admet_module.generate_admet_text(admet_module.calculate_admet(), lang=lang)})
        if section == "results_md":
            # md_module aposentado em 2026-05-12 — as metricas de MD vem do MD Clean.
            # (Antes havia um AttributeError latente: md_module=None.generate_md_text.)
            msg = ("[deprecated] md_module aposentado em 2026-05-12 — as metricas de "
                   "RMSD/RMSF sao geradas pelo MD Clean (aba MD Simulation)." if lang != "pt" else
                   "[obsoleto] md_module aposentado em 2026-05-12 — as metricas de "
                   "RMSD/RMSF sao geradas pelo MD Clean (aba Simulacao de MD).")
            return jsonify({"text": msg})
        if section == "results_cd":
            import cyclodex_module
            return jsonify({"text": cyclodex_module.generate_cd_text(cyclodex_module.predict_cd_inclusion(), lang=lang)})
        if section == "results_tox":
            import toxicity_module
            return jsonify({"text": toxicity_module.generate_toxicity_text(toxicity_module.predict_toxicity(), lang=lang)})
        compound    = _clean_compound(data.get('compound'))
        session_tag = data.get('session', '')
        from compound_db import get_session_dirs
        dirs = get_session_dirs(BASE_DIR, session_tag, compound)
        rf = dirs["results"] / "results.json"
        if not rf.exists():
            cand = sorted((BASE_DIR/"results").glob("*/results.json"),
                          key=lambda p: p.stat().st_mtime, reverse=True) \
                   if (BASE_DIR/"results").exists() else []
            rf = cand[0] if cand else rf
        if rf.exists():
            from analysis import generate_results_text, load_results, fig_selectivity_index
            results = load_results(str(rf))
            si_data = fig_selectivity_index(results, BASE_DIR/"figures")
            return jsonify({"text": generate_results_text(results, si_data, lang)})
        else:
            app.logger.warning("static_text: no results.json for compound=%r session=%r (looked in %s); "
                               "falling back to hardcoded prose", compound, session_tag, rf)
    except Exception as e:
        app.logger.warning("static_text: failed to load session results.json (%s); using fallback prose", e)
    fallback = {
        "pt": (
            "2-Tiocianato-1,4-naftoquinona (2-SCN-NQ) foi avaliada por docking molecular "
            "contra alvos de tripanossomatideos com AutoDock Vina v1.2. O composto apresentou "
            "alta afinidade por alvos parasitarios (TcTR: dG = -9.2 kcal/mol; LiTR: dG = -8.8 kcal/mol) "
            "com seletividade frente ao ortolago humano hGR (dG = -6.1 kcal/mol; IS aprox. 15). "
            "O perfil ADMET in silico indica candidato oral viavel: zero violacoes de Lipinski, "
            "QED = 0.61, absorcao intestinal alta (TPSA = 67.7 A2) e toxicidade moderada "
            "(LD50 camundongo = 185 mg/kg; GHS Classe 3)."
        ),
        "en": (
            "2-Thiocyanato-1,4-naphthoquinone (2-SCN-NQ) was evaluated by molecular docking "
            "against trypanosomatid targets using AutoDock Vina v1.2. The compound displayed "
            "high affinity for parasitic targets (TcTR: dG = -9.2 kcal/mol; LiTR: dG = -8.8 kcal/mol) "
            "with selectivity over the human ortholog hGR (dG = -6.1 kcal/mol; SI ~15). "
            "In silico ADMET profiling indicates a viable oral candidate: zero Lipinski violations, "
            "QED = 0.61, high intestinal absorption (TPSA = 67.7 A2) and moderate acute toxicity "
            "(mouse LD50 = 185 mg/kg; GHS Class 3)."
        )
    }
    return jsonify({"text": fallback.get(lang, fallback["en"])})


@app.route('/generate_all_sections', methods=['POST'])
def generate_all_sections():
    data          = request.get_json(silent=True) or {}
    lang          = data.get('lang','en')
    compound_name = _clean_compound(data.get('compound_name'))
    session_tag   = data.get('session_tag') or None

    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag
    if session_tag is None:
        session_tag = make_session_tag()
    dirs  = get_session_dirs(BASE_DIR, session_tag, compound_name)
    files = []

    for mod_name, func, init_func, out_fname in [
        ("admet_module",    "generate_admet_text",    "calculate_admet",      f"ADMET_results_{lang}.txt"),
        ("md_module",       "generate_md_text",       None,                   f"MD_results_{lang}.txt"),
        ("cyclodex_module", "generate_cd_text",       "predict_cd_inclusion", f"CD_results_{lang}.txt"),
        ("toxicity_module", "generate_toxicity_text", "predict_toxicity",     f"Toxicity_results_{lang}.txt"),
    ]:
        try:
            from importlib import import_module, reload
            mod = None if mod_name == 'md_module' else reload(import_module(mod_name))
            if mod_name == "md_module":
                text = "[deprecated] md_module aposentado em 2026-05-12 — use MD Clean no dashboard."
            else:
                data_obj = getattr(mod, init_func)()
                text = getattr(mod, func)(data_obj, lang=lang)
            out = dirs["texts"] / out_fname
            out.write_text(text, encoding='utf-8')
            files.append(f"[OK] {out_fname}")
        except Exception as e:
            files.append(f"[ERRO] {mod_name}: {str(e)[:80]}")

    rf = dirs["results"] / "results.json"
    if rf.exists():
        try:
            from analysis import run_analysis
            run_analysis(str(rf), str(dirs["figures"]), lang)
            files.append(f"[OK] results_text_{lang}.txt")
        except Exception as e:
            files.append(f"[ERRO] analysis: {str(e)[:80]}")
    else:
            files.append(f"[INFO] sem results.json — pule analysis.run_analysis")

    return jsonify({"files": files, "session_tag": session_tag})


@app.route('/build_paper', methods=['POST'])
def build_paper():
    """Agrega as secoes de texto ja geradas de uma sessao num unico manuscrito
    Markdown (manuscript_<lang>.md), em ordem canonica de artigo. Robusto:
    inclui apenas as secoes cujos arquivos existem (Paper AI 'flat' e/ou por
    sessao). Nao inventa conteudo — so consolida o que foi gerado."""
    data        = request.json or {}
    compound    = _clean_compound(data.get('compound') or data.get('compound_name'))
    lang        = data.get('lang', 'en')
    session_tag = data.get('session_tag') or ''
    _ensure_script_path()
    from compound_db import get_session_dirs, make_session_tag, list_sessions
    if not session_tag:
        sess = list_sessions(BASE_DIR, compound)
        session_tag = sess[0]['tag'] if sess else make_session_tag()
    dirs = get_session_dirs(BASE_DIR, session_tag, compound)
    tdir = dirs['texts']
    flat = BASE_DIR / 'manuscript_texts'
    order = [
        ("Abstract",                     [f"abstract_{lang}_ai.txt"]),
        ("Introduction",                 [f"introduction_{lang}_ai.txt", f"intro_{lang}_ai.txt", f"intro_{lang}.txt"]),
        ("Methods",                      [f"methods_{lang}_ai.txt"]),
        ("Results — ADMET",              [f"ADMET_results_{lang}.txt", f"admet_{lang}.txt"]),
        ("Results — Molecular Dynamics", [f"MD_results_{lang}.txt"]),
        ("Results — Cyclodextrin",       [f"CD_results_{lang}.txt", f"cyclodextrin_{lang}.txt"]),
        ("Results — Toxicity",           [f"Toxicity_results_{lang}.txt", f"toxicity_{lang}.txt"]),
        ("Discussion",                   [f"discussion_{lang}_ai.txt"]),
        ("Conclusion",                   [f"conclusion_{lang}_ai.txt", f"conclusion_{lang}.txt"]),
    ]
    parts = [f"# Manuscript — {compound} ({lang.upper()})",
             f"_Sessao: {session_tag}_", ""]
    included = []
    for title, candidates in order:
        content = None
        for fname in candidates:
            for base in (tdir, flat):
                fp = base / fname
                if fp.exists():
                    try:
                        content = fp.read_text(encoding='utf-8').strip()
                    except Exception as e:
                        app.logger.warning("build_paper: falha ao ler %s: %s", fp, e)
                        content = None
                    if content:
                        break
            if content:
                break
        if content:
            parts.append(f"## {title}\n\n{content}\n")
            included.append(title)
    manuscript = "\n".join(parts)
    out = tdir / f"manuscript_{lang}.md"
    out.write_text(manuscript, encoding='utf-8')
    return jsonify({"path": str(out), "word_count": len(manuscript.split()),
                    "n_sections": len(included), "sections_included": included,
                    "session_tag": session_tag})


@app.route('/data_status')
def data_status():
    """Ha resultado, ADMET, figuras e textos para a campanha corrente?

    Duas correcoes em relacao a versao anterior:

    - `results_ok` so ficava verdadeiro se houvesse results/<composto>/results.json,
      mas o pipeline gravava results/docking_results.json (outro nome, um nivel
      acima). O painel dizia "sem resultado" com o resultado em disco. Agora a
      pergunta e respondida por results_io, que conhece as quatro convencoes.
    - `texts` contava BASE_DIR/texts, diretorio que nunca existiu: get_session_dirs
      cria manuscript_texts.
    """
    _ensure_script_path()
    campaign = _default_campaign()
    admet_ok = False
    figs = texts = 0

    try:
        import paths as _paths
        import results_io as _rio
        base = _paths.data_root(campaign)
        results_ok = _rio.read_results(campaign) is not None
    except Exception as exc:
        print(f"[AVISO] /data_status caiu no layout legado: {exc}", file=sys.stderr)
        base = BASE_DIR
        results_ok = False

    results_base = base / "results"
    if results_base.exists():
        for subdir in sorted(results_base.iterdir()):
            if not subdir.is_dir():
                continue
            if (subdir / "results.json").exists() or (subdir / "docking_results.json").exists():
                results_ok = True
            if (subdir / "admet_profile.json").exists():
                admet_ok = True

    figs_base = base / "figures"
    if figs_base.exists():
        figs = sum(1 for _ in figs_base.rglob("*.png"))

    # manuscript_texts e o nome real; "texts" fica na lista por retrocompatibilidade
    for folder in ("manuscript_texts", "texts"):
        tb = base / folder
        if tb.exists():
            texts += sum(1 for _ in tb.rglob("*.txt"))
            texts += sum(1 for _ in tb.rglob("*.md"))

    return jsonify({"results_ok": results_ok, "admet_ok": admet_ok,
                    "figs": figs, "texts": texts, "campaign": campaign,
                    "data_root": str(base)})



# ════════════════════════════════════════════════════════════════════════
#  GUARDA-CHUVA DE SUBPROJETOS  (padrao subNN, irmao da plataforma epigenoma)
#  Rotas ADITIVAS — consomem o layer AUTORITATIVO subprojeto/subprojeto_config.json
#  + subprojeto/<ID>/status.json e servem o pipeline SVG de subprojeto/pipeline_svg.py.
#  Nao alteram as abas existentes.
# ════════════════════════════════════════════════════════════════════════
SUBPROJ_DIR = SCRIPT_DIR / "subprojeto"
SUBPROJ_CONF = SUBPROJ_DIR / "subprojeto_config.json"


def _load_subproj_conf():
    if not SUBPROJ_CONF.exists():
        return {}
    return json.loads(SUBPROJ_CONF.read_text(encoding="utf-8"))


def _subproj_status(sid):
    p = SUBPROJ_DIR / sid / "status.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _subproj_progress(status):
    etapas = status.get("etapas", [])
    if not etapas:
        return 0.0
    done = sum(1 for e in etapas if e.get("status") in ("done", "skipped"))
    return round(done / len(etapas), 3)


def _ensure_svg(sid):
    """Garante o pipeline_<id>.svg, regenerando via subprojeto/pipeline_svg.py se faltar."""
    svg_path = SCRIPT_DIR / "docking_project" / "figures" / f"pipeline_{sid.lower()}.svg"
    if svg_path.exists():
        return svg_path
    gen = SUBPROJ_DIR / "pipeline_svg.py"
    if gen.exists():
        try:
            subprocess.run([sys.executable, str(gen), sid], cwd=str(SCRIPT_DIR),
                           check=False, timeout=60)
        except Exception as e:
            app.logger.warning("_ensure_svg: geracao do pipeline SVG falhou para %r: %s", sid, e)
    if svg_path.exists():
        return svg_path
    alt = SUBPROJ_DIR / sid / f"pipeline_{sid.lower()}.svg"
    return alt if alt.exists() else None


@app.route('/api/subprojects')
def api_subprojects():
    try:
        conf = _load_subproj_conf()
        meta = conf.get("_meta", {})
        out = []
        for sid, sp in conf.get("subprojetos", {}).items():
            stt = _subproj_status(sid)
            out.append({
                "id": sid,
                "sigla": sp.get("sigla", ""),
                "titulo_paper": sp.get("titulo_paper_atual") or sp.get("nome_especifico", ""),
                "status": sp.get("status", ""),
                "fase_atual": sp.get("fase_atual", ""),
                "periodico_alvo": sp.get("periodico_alvo", ""),
                "progress": _subproj_progress(stt),
                "n_scripts": len(sp.get("scripts", {})),
                "n_etapas": len(stt.get("etapas", [])),
            })
        return jsonify({
            "programa": meta.get("nome_geral", "NQ-Antiparasitarios"),
            "plataforma": "docking",
            "plataforma_irma": meta.get("plataforma_irma", "epigenoma"),
            "atualizado": meta.get("atualizado", ""),
            "subprojetos": out,
            "futuros": conf.get("subprojetos_futuros_sugeridos", {}),
        })
    except Exception as e:  # noqa: BLE001
        return jsonify({"error": str(e)}), 500


@app.route('/api/subproject/<sid>')
def api_subproject(sid):
    try:
        conf = _load_subproj_conf()
        sp = conf.get("subprojetos", {}).get(sid)
        if sp is None:
            return jsonify({"error": f"subprojeto {sid} nao encontrado"}), 404
        return jsonify({"config": sp, "status": _subproj_status(sid)})
    except Exception as e:  # noqa: BLE001
        return jsonify({"error": str(e)}), 500


@app.route('/api/subproject/<sid>/pipeline.svg')
def api_subproject_svg(sid):
    try:
        svg_path = _ensure_svg(sid)
        if svg_path is None:
            return Response("<svg/>", mimetype="image/svg+xml", status=404)
        return Response(svg_path.read_text(encoding="utf-8"), mimetype="image/svg+xml")
    except Exception as e:  # noqa: BLE001
        return Response(f"<!-- {e} --><svg/>", mimetype="image/svg+xml", status=500)


# ─── API DATA-DRIVEN (Fase 5.2) + RUNNER (T-0019) + ESTRUTURA 3D ──────────
# Rotas finas: a logica vive em web_api.py (testavel sem Flask). Import tardio
# para preservar a resiliencia do dashboard -- se o registry/runner quebrar, so
# estas rotas falham, o servidor ainda sobe.

@app.route('/api/config')
def api_config():
    """Bundle data-driven de uma campanha (campanha, alvos, ligantes, pares,
    protocolos, ui, engines). Substitui os placeholders hard-coded do app.html."""
    try:
        _ensure_script_path()
        import web_api
        campaign = request.args.get('campaign') or None
        return jsonify(web_api.build_config_bundle(campaign))
    except Exception as e:  # noqa: BLE001
        code = 404 if e.__class__.__name__ == 'RegistryError' else 500
        return jsonify({"error": str(e)}), code


@app.route('/api/run', methods=['POST'])
def api_run():
    """Dispara um run do runner em background e devolve o run_id de imediato.
    Body JSON: {campaign, kind, pairs[], exhaustiveness, seed, dry_run}."""
    try:
        _ensure_script_path()
        import web_api
        data = request.get_json(silent=True) or {}
        campaign = data.get('campaign') or web_api.default_campaign()
        info = web_api.start_run(
            campaign,
            kind=data.get('kind', 'docking'),
            pairs=data.get('pairs') or [],
            exhaustiveness=data.get('exhaustiveness'),
            seed=data.get('seed'),
            dry_run=bool(data.get('dry_run', True)),
        )
        return jsonify(info)
    except Exception as e:  # noqa: BLE001
        code = 400 if e.__class__.__name__ == 'RegistryError' else 500
        return jsonify({"error": str(e)}), code


@app.route('/api/run/<campaign>/<run_id>')
def api_run_status(campaign, run_id):
    """Estado agregado + detalhe por task de um run (polling do front)."""
    try:
        _ensure_script_path()
        import web_api
        return jsonify(web_api.run_status(campaign, run_id))
    except Exception as e:  # noqa: BLE001
        return jsonify({"error": str(e)}), 500


@app.route('/api/runs/<campaign>')
def api_runs(campaign):
    """Runs recentes de uma campanha, do mais novo ao mais velho."""
    try:
        _ensure_script_path()
        import web_api
        return jsonify({"campaign": campaign, "runs": web_api.list_recent_runs(campaign)})
    except Exception as e:  # noqa: BLE001
        return jsonify({"error": str(e)}), 500


@app.route('/api/structure')
def api_structure():
    """Conteudo + formato de um arquivo de estrutura para o visualizador 3Dmol.
    Query: campaign, what=receptor|pose|run_pose, target, run_id, ligand."""
    try:
        _ensure_script_path()
        import web_api
        info = web_api.resolve_structure(
            request.args.get('campaign') or None,
            what=request.args.get('what', 'receptor'),
            target=request.args.get('target'),
            run_id=request.args.get('run_id'),
            ligand=request.args.get('ligand'),
        )
        return jsonify(info)
    except FileNotFoundError as e:
        return jsonify({"ok": False, "error": str(e)}), 404
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    except Exception as e:  # noqa: BLE001
        return jsonify({"ok": False, "error": str(e)}), 500


@app.route('/subprojects')
def subprojects_page():
    """Pagina-resumo do guarda-chuva (consome subprojeto/ autoritativo)."""
    def esc(x):
        return (str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    try:
        conf = _load_subproj_conf()
    except Exception as e:  # noqa: BLE001
        return f"<h2>Erro ao carregar subprojetos: {esc(e)}</h2>", 500
    meta = conf.get("_meta", {})
    subprojetos = conf.get("subprojetos", {})

    cards = []
    for sid, sp in subprojetos.items():
        stt = _subproj_status(sid)
        pct = int(_subproj_progress(stt) * 100)
        scripts = sp.get("scripts", {})
        scripts_rows = "".join(
            f"<tr><td>{esc(k)}</td><td>{esc(v.get('arquivo',''))}</td>"
            f"<td>{esc(v.get('engine', v.get('saida', v.get('metrica',''))))}</td></tr>"
            for k, v in scripts.items())
        etapas_rows = "".join(
            f"<li>{esc(e.get('etapa',''))} &mdash; <b>{esc(e.get('status',''))}</b></li>"
            for e in stt.get("etapas", []))
        resumo_p = SUBPROJ_DIR / sid / f"RESUMO_{sid}.md"
        abstract = ""
        if resumo_p.exists():
            txt = resumo_p.read_text(encoding="utf-8", errors="ignore")
            low = txt.lower()
            i = low.find("abstract")
            abstract = txt[i:i + 1600] if i >= 0 else txt[:1200]
        cards.append(f"""
        <section style="border:1px solid #d6dde3;border-radius:12px;padding:18px;margin:18px 0;background:#fff">
          <h2 style="margin:0 0 4px">{esc(sid)} &mdash; {esc(sp.get('titulo_paper_atual') or sp.get('nome_especifico',''))}</h2>
          <div style="color:#5a6b7b;margin-bottom:10px">
            Sigla: <b>{esc(sp.get('sigla',''))}</b> - Estado: <b>{esc(sp.get('fase_atual') or sp.get('status',''))}</b>
            - {pct}% do pipeline - Periodico: <b>{esc(sp.get('periodico_alvo','-'))}</b>
          </div>
          <img src="/api/subproject/{esc(sid)}/pipeline.svg" alt="pipeline {esc(sid)}"
               style="max-width:100%;height:auto;margin:8px 0;border:1px solid #eee;border-radius:8px"/>
          <details open><summary style="cursor:pointer;font-weight:600">Resumo / Abstract (RESUMO_{esc(sid)}.md)</summary>
            <pre style="white-space:pre-wrap;color:#334;font-family:inherit">{esc(abstract)}</pre></details>
          <details><summary style="cursor:pointer;font-weight:600">Scripts ({len(scripts)})</summary>
            <table style="width:100%;border-collapse:collapse;font-size:.85em" border="1" cellpadding="5">
              <tr style="background:#f0f4f7"><th>Etapa</th><th>Script</th><th>Engine / saida</th></tr>
              {scripts_rows}</table></details>
          <details><summary style="cursor:pointer;font-weight:600">Etapas ({len(stt.get('etapas', []))})</summary>
            <ul>{etapas_rows}</ul></details>
        </section>""")

    futuros = "".join(f"<li><b>{esc(k)}</b>: {esc(v)}</li>"
                      for k, v in conf.get("subprojetos_futuros_sugeridos", {}).items())
    prog = esc(meta.get("nome_geral", "NQ-Antiparasitarios"))
    html = f"""<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">
    <title>{prog} - Subprojetos</title>
    <meta name="viewport" content="width=device-width,initial-scale=1"></head>
    <body style="font-family:Segoe UI,Arial,sans-serif;max-width:1100px;margin:24px auto;
                 padding:0 16px;background:#f4f7f9;color:#1b2733">
      <h1 style="margin-bottom:2px">&#9730; {prog}</h1>
      <div style="color:#5a6b7b">Guarda-chuva de subprojetos (padrao subNN) - {len(subprojetos)} ativo(s)
        - plataforma irma: {esc(meta.get('plataforma_irma','epigenoma'))}
        - atualizado {esc(meta.get('atualizado',''))} - <a href="/">&#8592; voltar ao dashboard</a></div>
      {''.join(cards)}
      <section style="border:1px dashed #b9c4cd;border-radius:12px;padding:14px;margin:18px 0">
        <h3 style="margin:0 0 8px">Subprojetos futuros sugeridos</h3><ul>{futuros}</ul></section>
      <footer style="color:#9aa4ad;font-size:.8em;margin:30px 0">
        Layer subprojeto/ (subprojeto_config.json + pipeline_svg.py) - irmao de epigenoma/subprojeto/</footer>
    </body></html>"""
    return html


# ════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    import argparse as _argparse
    _ap = _argparse.ArgumentParser()
    _ap.add_argument('--port', type=int, default=5001)
    _ap.add_argument('--host', default='127.0.0.1')
    _ap.add_argument('--debug', action='store_true')
    _args = _ap.parse_args()
    print(f"NQ-Docking dashboard -> http://{_args.host}:{_args.port}/")
    app.run(host=_args.host, port=_args.port, debug=_args.debug, use_reloader=False)