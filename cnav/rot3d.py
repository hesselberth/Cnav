#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 17:23:15 2023

@author: Marcel Hesselberth
"""

# =============================================================================
# 3D rotations following the sign convention in the astronomical almanac.
# =============================================================================


from numpy import array, sin, cos, double


def R1(phi, floattype=double):
    s = sin(phi)
    c = cos(phi)
    return array([[1, 0, 0], [0, c, s], [0, -s, c]], dtype=floattype)

def R2(phi, floattype=double):
    s = sin(phi)
    c = cos(phi)
    return array([[c, 0, -s], [0, 1, 0], [s, 0, c]], dtype=floattype)

def R3(phi, floattype=double):
    s = sin(phi)
    c = cos(phi)
    return array([[c, s, 0], [-s, c, 0], [0, 0, 1]], dtype=floattype)
