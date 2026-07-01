#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""conftest.py - garante que o pacote 'bonder' esteja importavel nos testes."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
