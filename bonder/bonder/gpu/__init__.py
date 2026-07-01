#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bonder.gpu - aceleracao GPU (RTX 3060) com fallback CPU garantido.

detect.probe()  : escolha segura de device (cuda/cpu) com teste de kernel real.
gpu_gate        : lease cooperativo com o LLM local (Ollama).
tools           : DIA-NN / MS2Rescore / ColabFold (consagrados, revisados).
esm2            : priorizacao exploratoria (SUPLEMENTAR, off por padrao).
"""
from bonder.gpu.detect import probe, DeviceDecision, reset_cache
from bonder.gpu import gpu_gate, tools

__all__ = ["probe", "DeviceDecision", "reset_cache", "gpu_gate", "tools"]
