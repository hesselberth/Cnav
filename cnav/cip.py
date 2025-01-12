#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 20:28:19 2024

@author: Marcel Hesselberth
"""

import numpy as np
from .constants import AS2RAD, UAS2RAD, PI, PI2, DEG2RAD
from .rot3d import *
from .xys import XYs06 as XYs, PFW06_gamma_phi as PFW


# Celestial to terrestrial coordinate transformations based on the X, Y
# components of the CIP unit vector and on quantity s, the CIO locator.

# ToDo: check units. All angles radians.


def Mcio(tjc):
    M = np.empty((3,3))

    X, Y, s = XYs(tjc)
    Z = np.sqrt(1 - X*X - Y*Y)
    a = 1 / (1+Z)
    ss, cs = np.sin(s), np.cos(s)

    M[0][0] = cs + a * X * ( Y * ss - X * cs)
    M[0][1] = -ss + a * Y * (Y * ss - X * cs)
    M[0][2] = -(X * cs - Y * ss)
    M[1][0] = ss - a * X * (Y * cs + X * ss)
    M[1][1] = cs - a * Y * (Y * cs + X * ss)
    M[1][2] = -(Y * cs + X * ss)
    M[2][0] = X
    M[2][1] = Y
    M[2][2] = Z
    return M


def ERA(UT1_2000):  # UT1 is Julian UT1 date since JD2000
    T_u = UT1_2000
    f = UT1_2000 % 1  # TODO: check for negative UT1_2000
    turns = (f + 0.7790572732640 + 0.00273781191135448 * T_u) % 1
    return PI2 * turns


def R(ut1_2000, tjc, mcio = None):
    if not mcio:
        mcio = Mcio(tjc)
    return R3(ERA(ut1_2000)) @ mcio


# classical equinox based NPB matrix and EO
def Mclass_EO(tjc):
    gamma, phi = PFW(tjc)
    k = np.array([np.sin(phi)*np.sin(gamma), \
                  -np.sin(phi)*np.cos(gamma), \
                  np.cos(phi)])
    X, Y, s = XYs(tjc)
    Z = np.sqrt(1 - X*X - Y*Y)
    n = np.array([X, Y, Z])
    nxk = np.cross(n, k)
    nxnxk = np.cross(n, nxk)
    YY = nxk / np.linalg.norm(nxk)
    y = nxnxk / np.linalg.norm(nxnxk)
    Mclass = np.array([YY, y, n], dtype=np.double)
    
    zp1 = 1 + Z
    S = np.array([1 - X * X / zp1, -X * Y / zp1, -X], dtype = np.double)

    return Mclass, s - np.arctan2(np.dot(y, S), np.dot(YY, S))


def Mclass(Mcio, tjc):
    gamma, phi = PFW(tjc)
    k = np.array([np.sin(phi)*np.sin(gamma), \
                  -np.sin(phi)*np.cos(gamma), \
                  np.cos(phi)])
    n = Mcio[2]
    nxk = np.cross(n, k)
    nxnxk = np.cross(n, nxk)
    YY = nxk / np.linalg.norm(nxk)
    y = nxnxk / np.linalg.norm(nxnxk)
    Mclass = np.array([YY, y, n], dtype=np.double)
    return Mclass
 
    
def EO(tjc):
    X, Y, S = XYs(tjc)
    return Mclass_EO(tjc)[1]


def GST(UT1_2000, eo):
        return ERA(UT1_2000) - eo

