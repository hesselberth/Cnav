#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 05 16:27:27 2025

@author: Marcel Hesselberth
"""

import os
import sys
sys.path.append(os.getcwd() + '/..')


def test_numba_installed():
    import numba

def test_cnumba():
    import cnav.cnumba

def test_acc():
    import cnav.cnumba
    assert(cnav.cnumba.numba_acc)

