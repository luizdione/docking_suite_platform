# -*- coding: utf-8 -*-
"""
engines/boltz2.py -- motor Boltz-2 (predicao de estrutura + afinidade).

Boltz-2 preve o complexo proteina-ligante e uma estimativa de afinidade. Precisa
de um YAML de entrada com a SEQUENCIA da proteina (que o registry nao guarda -- e
metadado estrutural, nao de sequencia) e o ligante. A sequencia vem de
req.inputs['sequence'] ou de um FASTA em req.inputs['fasta']; o SMILES do ligante
sai do registry. Sem sequencia -> erro claro (nunca inventa).

O YAML e escrito a mao (sem depender de PyYAML). A execucao real dispara
`boltz predict` (o binario foi resolvido no env 'foldgpu'); dry_run so planeja.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from engines import toolchain
from engines.base import Engine, EngineRequest, EngineResult, register


def _ligand_smiles(req):
    smiles = (req.params or {}).get("smiles")
    if smiles:
        return smiles
    try:
        import registry
        return registry.get_ligand(req.ligand).get("smiles")
    except Exception:
        return None


def _sequence(req):
    seq = (req.inputs or {}).get("sequence")
    if seq:
        return "".join(str(seq).split())
    fasta = (req.inputs or {}).get("fasta")
    if fasta and Path(fasta).exists():
        lines = Path(fasta).read_text(encoding="utf-8").splitlines()
        return "".join(ln.strip() for ln in lines if ln and not ln.startswith(">"))
    return None


def _yaml(sequence, smiles, want_affinity):
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {sequence}",
             "  - ligand:", "      id: L", f"      smiles: {smiles}"]
    if want_affinity:
        lines += ["properties:", "  - affinity:", "      binder: L"]
    return "\n".join(lines) + "\n"


@register
class Boltz2Engine(Engine):
    name = "boltz2"

    def capabilities(self):
        return {"pose", "affinity", "structure", "confidence"}

    def check_toolchain(self):
        boltz = toolchain.resolve_known("boltz")
        return {"available": bool(boltz),
                "tools": {"boltz": toolchain.tool_entry(boltz) if boltz
                          else {"available": False}},
                "gpu": toolchain.gpu_entry()}

    def prepare(self, req: EngineRequest) -> EngineRequest:
        seq = _sequence(req)
        smiles = _ligand_smiles(req)
        if not seq or not smiles:
            return req  # run() reporta o erro com detalhe
        content = _yaml(seq, smiles, "affinity" in self.capabilities())
        if req.workdir and not req.dry_run:
            Path(req.workdir).mkdir(parents=True, exist_ok=True)
            yaml_path = Path(req.workdir) / f"boltz_{req.target}__{req.ligand}.yaml"
            yaml_path.write_text(content, encoding="utf-8")
            req.inputs["yaml"] = str(yaml_path)
        req.inputs["_yaml_preview"] = content
        return req

    def _cmdline(self, req):
        p = req.protocol.get("params") if req.protocol else {}
        p = p or {}
        cmd = ["boltz", "predict", req.inputs.get("yaml", "<input.yaml>")]
        if req.workdir:
            cmd += ["--out_dir", str(req.workdir)]
        cmd += ["--diffusion_samples", str(p.get("diffusion_samples", 5))]
        cmd += ["--output_format", str(p.get("output_format", "mmcif"))]
        if p.get("use_msa_server"):
            cmd.append("--use_msa_server")
        if p.get("write_full_pae"):
            cmd.append("--write_full_pae")
        return cmd

    def run(self, req: EngineRequest) -> EngineResult:
        req = self.prepare(req)
        seq = _sequence(req)
        smiles = _ligand_smiles(req)
        cmd = self._cmdline(req)
        if not seq:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error="sequencia da proteina ausente: forneca "
                                      "inputs['sequence'] ou inputs['fasta'] "
                                      "(o registry nao guarda sequencia)")
        if not smiles:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"SMILES do ligante '{req.ligand}' nao resolvido")
        if req.dry_run:
            return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                                metrics={"mode": "dry_run",
                                         "yaml": req.inputs.get("_yaml_preview")})
        boltz = toolchain.resolve_known("boltz")
        if not boltz:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error="boltz nao encontrado (env conda 'foltgpu'/'foldgpu' "
                                      "com boltz instalado)")
        cmd[0] = boltz
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=req.params.get("timeout_s", 86400))
        except Exception as exc:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"boltz subprocess falhou: {exc!r}")
        if proc.returncode != 0:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"boltz exit {proc.returncode}",
                                log=(proc.stderr or proc.stdout or "")[-800:])
        return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                            outputs={"out_dir": str(req.workdir) if req.workdir else None},
                            metrics={"returncode": 0}, log=(proc.stdout or "")[-2000:])
