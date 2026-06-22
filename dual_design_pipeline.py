# -*- coding: utf-8 -*-
"""
dual_design_pipeline.py  —  Engenharia DUAL (substrato subversivo + inibidor classico)
da 2-tiocianil-1,4-naftoquinona contra a TRIPANOTIONA REDUTASE (TcTR / Leishmania TR).

Objetivo: propor 5 estruturas MAIORES, derivadas da naftoquinona, que exercam
SIMULTANEAMENTE as duas funcoes:
  (A) SUBSTRATO SUBVERSIVO — o nucleo naftoquinonico (quinona) compete por eletrons
      do FAD/FADH-, e ciclo-reduz gerando superoxido/radicais livres (estresse oxidativo).
  (B) INIBIDOR CLASSICO — uma CAUDA acoplada ocupa o canal/entrada do substrato
      (sitio hidrofobico/cationico tipo "mepacrina"; o substrato natural tripanotiona
      tem espinha de espermidina), bloqueando a ligacao da tripanotiona (competitivo).

Roda LOCAL (RDKit + opcional ollama), em background, SEM custo de tokens externos.
Saidas (em docking_project/figures/dual_design/):
  - dual_design_results.json        (dados estruturados)
  - dual_design_report.html         (relatorio self-contained com imagens 2D)
  - struct_<id>.png                 (estruturas)
Uso: python dual_design_pipeline.py [--ollama]
"""
import os, sys, json, base64, datetime, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "docking_project", "figures", "dual_design")
os.makedirs(OUTDIR, exist_ok=True)
LOG = os.path.join(OUTDIR, "dual_design.log")


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


PARENT = {"name": "2-SCN-NQ (parental)", "smiles": "O=C1C=C(SC#N)C(=O)c2ccccc21"}

# ── ESTADO DA ARTE (TcTR / Leishmania TR) — sintese curada (medchem real) ────────
STATE_OF_ART = """\
TR como alvo (TcTR e Leishmania TR/LiTR): enzima homodimerica FAD-dependente, exclusiva
de tripanossomatideos (humano usa glutationa redutase, hGR). Dois sub-sitios exploraveis:
(1) o sitio CATALITICO/FAD com o ditiol redox Cys52/Cys57 — onde substratos subversivos
(quinonas, nitrofuranos, fenotiazinas) sao reduzidos pelo FADH- e ciclam gerando ERO;
(2) um amplo CANAL hidrofobico/anionico de entrada do substrato (tripanotiona = bis-
glutationil-espermidina), que acomoda POLIAMINAS e cations lipofilicos — o "sitio
mepacrina/Z" onde se ligam inibidores classicos competitivos.
Estrategias consagradas: (a) conjugados POLIAMINA-quinona/aciridina (a poliamina mimetiza
a espermidina e ancora no canal; ex. revisoes de Krauth-Siegel; Cavalli & Bolognesi 2009);
(b) FENOTIAZINAS (tioridazina, clorpromazina) e ACRIDINAS/mepacrina como inibidores do
sitio hidrofobico; (c) serie NQ-TRIAZOL (Bahia 2016) e NQ-2-substituidas (De Castro/
Menna-Barreto) como redox-ciclantes; (d) HIBRIDOS/bivalentes que cruzam os dois sub-sitios.
Seletividade vs hGR vem da arquitetura do canal (TcTR mais largo/menos cationico; Ala34
onde hGR tem Arg) — favorece ligantes volumosos e poliamina-mimeticos.
LACUNA que este projeto explora: o farmaco atual (2-SCN-NQ) liga ESTAVEL no sitio de
ENTRADA (~9 A do FAD; modo INIBIDOR, confirmado por docking+MD), mas e PEQUENO demais
para tocar o FAD. AUMENTAR a molecula com um nucleo redox posicionado para o FAD + cauda
de canal => molecula DUAL (subversivo no FAD + inibidor competitivo no canal)."""

# ── 5 DESENHOS DUAIS (curados, fundamentados) ────────────────────────────────────
# Cada um: nucleo NQ (redox/subversivo) + linker + grupo de canal (inibidor classico),
# MAIOR que o parental, projetado para alcancar AMBOS os sub-sitios.
DESIGNS = [
    {
        "id": "NQ-SPM",
        "name": "NQ-2-tio-espermidina (conjugado poliamina)",
        "smiles": "O=C1C=C(SCCCNCCCCNCCCN)C(=O)c2ccccc21",
        "dual": "Nucleo 1,4-NQ = redox/subversivo no FAD; cauda espermidina (3-4-3 triamina) "
                "mimetiza a tripanotiona e ancora no canal anionico de entrada (inibidor competitivo).",
        "why_bigger": "Estende ~8-12 A em direcao ao canal, alcancando ambos os sub-sitios; "
                      "poliamina protonada da seletividade vs hGR (canal mais estreito/cationico).",
        "synthesis": "1) 2-bromo-1,4-NQ + N1,N4-bis(Boc)-N8-(3-mercaptopropil)espermidina, "
                     "SNAr tiol-NQ (K2CO3/DMF, 60 C); 2) desprotecao Boc (TFA/DCM). "
                     "Alt.: tiol-eno do 2-vinil-NQ com tiol-espermidina.",
        "refs": "Conjugados poliamina-quinona vs TR (Krauth-Siegel; Cavalli & Bolognesi 2009).",
    },
    {
        "id": "NQ-PTZ",
        "name": "NQ-fenotiazina (hibrido redox + inibidor de sitio hidrofobico)",
        "smiles": "O=C1C=C(NCCCN2c3ccccc3Sc3ccccc32)C(=O)c2ccccc21",
        "dual": "NQ redox no FAD; fenotiazina (farmacoforo de TR validado: tioridazina/"
                "clorpromazina) ocupa o bolso hidrofobico de entrada (inibidor classico).",
        "why_bigger": "Fenotiazina volumosa preenche o sitio mepacrina; tambem mitocondriotropica "
                      "(catio lipofilico) — reforca o estresse oxidativo no parasito.",
        "synthesis": "10-(3-aminopropil)fenotiazina + 2-bromo-1,4-NQ, SNAr (K2CO3/DMF, 80 C, 6-12 h). "
                     "Aminopropil-fenotiazina vem de fenotiazina + 1-bromo-3-cloropropano/NaH e aminacao.",
        "refs": "Fenotiazinas como inibidores de TR (Krauth-Siegel; Khan 2007).",
    },
    {
        "id": "NQ-ACR",
        "name": "NQ-acridina (mepacrina-mimetico)",
        "smiles": "O=C1C=C(NCCCNc3c4ccccc4nc4ccccc34)C(=O)c2ccccc21",
        "dual": "NQ redox no FAD; 9-aminoacridina (mepacrina-mimetico) ancora no sitio "
                "hidrofobico/Z de entrada da tripanotiona — inibidor competitivo classico.",
        "why_bigger": "Acridina planar grande intercala/empacota no canal; linker propil ajusta "
                      "o alcance ate o FAD para o nucleo NQ.",
        "synthesis": "1) 9-cloroacridina + N-Boc-1,3-propanodiamina (fenol fundido, 100 C) -> "
                     "9-(3-aminopropilamino)acridina; 2) SNAr com 2-bromo-1,4-NQ; 3) (sem Boc aqui).",
        "refs": "Mepacrina/acridinas como inibidores de TR (Bonse 1999; Saravanamuthu 2004).",
    },
    {
        "id": "NQ-TZ-POL",
        "name": "NQ-triazol-poliamina (click, modular)",
        "smiles": "O=C1C=C(n2cc(CNCCCNCCCCN)nn2)C(=O)c2ccccc21",
        "dual": "Aproveita a serie NQ-triazol (redox, Bahia 2016); o braco poliamina via "
                "triazol-CH2 alcanca o canal de entrada (inibidor). Triazol adiciona H-bond/pi perto do FAD.",
        "why_bigger": "Modularidade do CuAAC permite ajustar o comprimento da poliamina para "
                      "vencer a distancia FAD<->canal (~9-12 A).",
        "synthesis": "CuAAC: 2-azido-1,4-NQ + N-propargil-poliamina(Boc) (CuSO4/Na-ascorbato, "
                     "t-BuOH/H2O); depois desprotecao Boc (TFA).",
        "refs": "Serie NQ-triazol (Bahia et al. 2016 MedChemComm 7:1555).",
    },
    {
        "id": "BIS-NQ",
        "name": "Bis-naftoquinona (dois nucleos redox, ponte de canal)",
        "smiles": "O=C1C=C(SCCNCCCCNCCSC3=CC(=O)c4ccccc4C3=O)C(=O)c2ccccc21",
        "dual": "Dois nucleos NQ ligados por diamina/ditioeter: um voltado ao FAD (subversivo, "
                "dobra a capacidade redox) e outro ocupando a entrada; a diamina central ancora no canal.",
        "why_bigger": "Bivalente: cruza os dois sub-sitios de uma vez; aumenta avidez e o fluxo "
                      "de eletrons (mais ERO).",
        "synthesis": "2 equiv 2-bromo-1,4-NQ + N,N'-bis(2-mercaptoetil)-1,4-butanodiamina, "
                     "dupla SNAr tiol-NQ (K2CO3/DMF). Estequiometria controla mono vs bis.",
        "refs": "Bis-quinonas/ligantes bivalentes de TR; conceito de inibidor dual-site.",
    },
]


# ── Docking automatico nos DOIS sub-sitios (Vina local, CPU — nao compete c/ a GPU do MD) ──
VINA = os.path.join(HERE, "vina.exe")
RECEPTOR = os.path.join(HERE, "docking_project", "proteins", "1BZL_receptor.pdbqt")
SITES = {
    # sitio CATALITICO/redox (ditiol Cys52/57 a ~4 A do FAD) -> engajamento SUBVERSIVO
    "FAD": dict(cx=28.02, cy=6.75, cz=-11.34, sx=24, sy=24, sz=24),
    # canal de ENTRADA do substrato (onde o parental 2-SCN-NQ liga) -> INIBIDOR competitivo
    "ENTRADA": dict(cx=47.4, cy=-0.5, cz=3.7, sx=22, sy=22, sz=22),
}


OBABEL = os.path.join(r"C:\Users\luizd\miniforge3\envs\datascience\Library\bin", "obabel.exe")


def smiles_to_pdbqt(smiles, out_pdbqt):
    """SMILES -> 3D via RDKit (ETKDG+MMFF, alta qualidade) -> pdbqt via obabel CLI
    (so converte; sem 3D-gen do obabel, evitando o problema do UFF.prm)."""
    import subprocess
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return False
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
            AllChem.EmbedMolecule(m, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:
            pass
        pdb = out_pdbqt.replace(".pdbqt", "_3d.pdb")
        Chem.MolToPDBFile(m, pdb)
        obabel = OBABEL if os.path.exists(OBABEL) else "obabel"
        subprocess.run([obabel, pdb, "-opdbqt", "-O", out_pdbqt],
                       capture_output=True, text=True, timeout=120)
        return os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0
    except Exception as e:
        log("  (prep pdbqt falhou p/ %s: %s)" % (smiles[:30], e))
        return False


def dock_site(lig_pdbqt, site, out_pdbqt):
    """Roda Vina (exhaustiveness 16, 4 CPUs) e devolve a melhor afinidade (kcal/mol)."""
    import subprocess
    try:
        subprocess.run([VINA, "--receptor", RECEPTOR, "--ligand", lig_pdbqt,
                        "--center_x", str(site["cx"]), "--center_y", str(site["cy"]),
                        "--center_z", str(site["cz"]), "--size_x", str(site["sx"]),
                        "--size_y", str(site["sy"]), "--size_z", str(site["sz"]),
                        "--exhaustiveness", "16", "--cpu", "4", "--seed", "42",
                        "--out", out_pdbqt],
                       capture_output=True, text=True, timeout=1200)
        if os.path.exists(out_pdbqt):
            for ln in open(out_pdbqt, errors="ignore"):
                if ln.startswith("REMARK VINA RESULT"):
                    return float(ln.split(":")[1].split()[0])
    except Exception as e:
        log("  (dock falhou: %s)" % e)
    return None


def main():
    use_ollama = "--ollama" in sys.argv
    do_dock = ("--no-dock" not in sys.argv) and os.path.exists(VINA) and os.path.exists(RECEPTOR)
    open(LOG, "w").close()
    log("=== Dual-design pipeline (subversivo + inibidor) — TcTR/Leishmania TR ===")
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Draw, rdMolDescriptors, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception as e:
        log("ERRO: RDKit indisponivel (%s). Rode no env datascience." % e)
        sys.exit(2)

    def props(smiles):
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None, None
        d = dict(
            mw=round(Descriptors.MolWt(m), 1),
            clogp=round(Descriptors.MolLogP(m), 2),
            tpsa=round(Descriptors.TPSA(m), 1),
            hbd=rdMolDescriptors.CalcNumHBD(m),
            hba=rdMolDescriptors.CalcNumHBA(m),
            rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
            arom=rdMolDescriptors.CalcNumAromaticRings(m),
            nbasic=sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0),
        )
        return m, d

    def draw(m, path):
        try:
            AllChem.Compute2DCoords(m)
            dr = rdMolDraw2D.MolDraw2DCairo(420, 320)
            dr.DrawMolecule(m)
            dr.FinishDrawing()
            open(path, "wb").write(dr.GetDrawingText())
            return True
        except Exception as e:
            log("  (draw falhou: %s)" % e)
            return False

    pm, pprops = props(PARENT["smiles"])
    log("parental %s: MW=%s clogp=%s" % (PARENT["name"], pprops["mw"], pprops["clogp"]))
    log("docking automatico (Vina, 2 sub-sitios): %s" % ("ON" if do_dock else "OFF"))

    results = []
    for d in DESIGNS:
        m, pr = props(d["smiles"])
        if m is None:
            log("  SMILES invalido em %s — pulando." % d["id"])
            continue
        png = os.path.join(OUTDIR, "struct_%s.png" % d["id"])
        ok = draw(m, png)
        img_b64 = base64.b64encode(open(png, "rb").read()).decode() if ok else ""
        grew = pr["mw"] > pprops["mw"] + 80
        has_quinone = m.HasSubstructMatch(Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")) or \
                      m.HasSubstructMatch(Chem.MolFromSmarts("O=C1C=CC(=O)C=C1"))
        has_channel = pr["nbasic"] >= 1 or pr["arom"] >= 3

        # ── DOCKING nos dois sub-sitios ──
        gFAD = gENT = dual_dG = None
        if do_dock:
            ligp = os.path.join(OUTDIR, "lig_%s.pdbqt" % d["id"])
            if smiles_to_pdbqt(d["smiles"], ligp):
                gFAD = dock_site(ligp, SITES["FAD"], os.path.join(OUTDIR, "dock_%s_FAD.pdbqt" % d["id"]))
                gENT = dock_site(ligp, SITES["ENTRADA"], os.path.join(OUTDIR, "dock_%s_ENT.pdbqt" % d["id"]))
                log("    %s dock: FAD=%s  ENTRADA=%s" % (d["id"], gFAD, gENT))

        # ── SCORE: dirigido por afinidade quando disponivel ──
        # Um bom ligante DUAL liga BEM nos DOIS sub-sitios; rank por -media(dG),
        # bonus se ambos fortes (< -7) e proximos (atinge os dois de fato).
        druglike_pen = (pr["mw"] > 700) + (pr["clogp"] > 6 or pr["clogp"] < -1)
        if gFAD is not None and gENT is not None:
            dual_dG = round((gFAD + gENT) / 2.0, 2)
            both_strong = (gFAD < -7.0 and gENT < -7.0)
            balanced = abs(gFAD - gENT) < 1.5
            score = round(-dual_dG + 1.0 * both_strong + 0.5 * balanced
                          + 1.0 * bool(has_quinone) + 1.0 * has_channel - druglike_pen, 2)
        else:
            # fallback heuristico (sem docking)
            score = round(2.0 * grew + 2.0 * bool(has_quinone) + 2.0 * has_channel - druglike_pen, 2)

        d2 = dict(d)
        d2.update(pr)
        d2["dMW_vs_parent"] = round(pr["mw"] - pprops["mw"], 1)
        d2["redox_core_ok"] = bool(has_quinone)
        d2["channel_anchor_ok"] = bool(has_channel)
        d2["dG_FAD"] = gFAD
        d2["dG_ENTRADA"] = gENT
        d2["dual_dG"] = dual_dG
        d2["dual_score"] = score
        d2["img"] = img_b64
        results.append(d2)
        log("  %s: MW=%s (+%s) clogp=%s quinona=%s canal=%s dualdG=%s score=%s"
            % (d["id"], pr["mw"], d2["dMW_vs_parent"], pr["clogp"], has_quinone,
               has_channel, dual_dG, score))

    results.sort(key=lambda r: -r["dual_score"])

    narrative = ""
    if use_ollama:
        narrative = _ollama_narrative(results)

    payload = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "parent": {**PARENT, **pprops},
        "state_of_art": STATE_OF_ART,
        "designs": results,
        "ollama_narrative": narrative,
        "n": len(results),
    }
    json_path = os.path.join(OUTDIR, "dual_design_results.json")
    json.dump(payload, open(json_path, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log("JSON -> %s" % json_path)

    _write_html(payload, os.path.join(OUTDIR, "dual_design_report.html"))
    log("HTML -> %s" % os.path.join(OUTDIR, "dual_design_report.html"))
    log("=== concluido: %d estruturas duais ===" % len(results))


def _ollama_narrative(results):
    """Opcional: usa o LLM local (ollama) p/ ampliar o racional. Silencioso se indisponivel."""
    try:
        import urllib.request
        names = "; ".join("%s (%s)" % (r["name"], r["smiles"]) for r in results)
        prompt = ("Voce e um quimico medicinal. Em 200 palavras, em portugues, comente a "
                  "estrategia de tornar a 2-tiocianil-1,4-naftoquinona um ligante DUAL "
                  "(substrato subversivo no FAD + inibidor competitivo no canal da tripanotiona) "
                  "da tripanotiona redutase, dados estes candidatos: " + names)
        req = urllib.request.Request(
            "http://localhost:11434/api/generate",
            data=json.dumps({"model": "llama3", "prompt": prompt, "stream": False}).encode(),
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode()).get("response", "").strip()
    except Exception as e:
        log("  (ollama narrativa pulada: %s)" % e)
        return ""


def _write_html(payload, path):
    P = payload["parent"]
    rows = []
    for i, d in enumerate(payload["designs"], 1):
        img = ('<img src="data:image/png;base64,%s" style="width:340px;border:1px solid #ccc;border-radius:6px">' % d["img"]) if d.get("img") else ""
        rows.append("""
        <div style="border:1px solid #d0d7de;border-radius:10px;padding:14px;margin:12px 0;background:#fff">
          <h3 style="margin:0 0 6px">#%d &middot; %s <span style="color:#0969da">[score %.1f]</span></h3>
          <div style="display:flex;gap:16px;flex-wrap:wrap">
            <div>%s</div>
            <div style="flex:1;min-width:280px;font-size:.9em;line-height:1.5">
              <code style="background:#f6f8fa;padding:2px 5px;border-radius:4px">%s</code><br>
              <b>MW</b> %.1f (Δ +%.1f vs parental) &middot; <b>cLogP</b> %.2f &middot; <b>TPSA</b> %.1f &middot;
              <b>HBD/HBA</b> %d/%d &middot; <b>rot</b> %d &middot; <b>aneis arom.</b> %d<br>
              <b>Docking ΔG (kcal/mol):</b> FAD/redox %s &middot; ENTRADA/canal %s &middot; média(dual) %s<br>
              <b>Nucleo redox (subversivo):</b> %s &middot; <b>Ancora de canal (inibidor):</b> %s<br>
              <b>Funcao dual:</b> %s<br>
              <b>Por que maior:</b> %s<br>
              <b>Sintese sugerida:</b> %s<br>
              <span style="color:#57606a"><b>Lit.:</b> %s</span>
            </div>
          </div>
        </div>""" % (i, d["name"], d["dual_score"], img, d["smiles"], d["mw"], d["dMW_vs_parent"],
                     d["clogp"], d["tpsa"], d["hbd"], d["hba"], d["rotb"], d["arom"],
                     ("%.2f" % d["dG_FAD"]) if d.get("dG_FAD") is not None else "n/d",
                     ("%.2f" % d["dG_ENTRADA"]) if d.get("dG_ENTRADA") is not None else "n/d",
                     ("%.2f" % d["dual_dG"]) if d.get("dual_dG") is not None else "n/d",
                     "OK" if d["redox_core_ok"] else "-", "OK" if d["channel_anchor_ok"] else "-",
                     d["dual"], d["why_bigger"], d["synthesis"], d["refs"]))
    narr = ("<h2>Racional (LLM local)</h2><p style='white-space:pre-wrap'>%s</p>" % payload["ollama_narrative"]) if payload["ollama_narrative"] else ""
    html = """<!doctype html><html lang="pt"><head><meta charset="utf-8">
<title>Engenharia Dual — NQ vs Tripanotiona Redutase</title></head>
<body style="font-family:Inter,Arial,sans-serif;max-width:980px;margin:18px auto;color:#1f2328">
<h1>&#129516; Engenharia DUAL da naftoquinona &mdash; subversivo + inibidor (TcTR / Leishmania TR)</h1>
<p style="color:#57606a">Gerado %s &middot; parental %s (MW %.1f). 5 estruturas maiores projetadas
para tocar o FAD (redox/subversivo) E ocupar o canal de entrada da tripanotiona (inibidor classico).</p>
<div style="background:#f6f8fa;border-left:4px solid #0969da;padding:10px 14px;font-size:.9em;white-space:pre-wrap;line-height:1.5">%s</div>
%s
<h2>5 estruturas modificadas (ordenadas por score dual)</h2>
%s
</body></html>""" % (payload["generated"], P["name"], P["mw"], payload["state_of_art"], narr, "".join(rows))
    open(path, "w", encoding="utf-8").write(html)


if __name__ == "__main__":
    main()
