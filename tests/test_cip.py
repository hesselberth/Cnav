#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 18:38:37 2024

@author: Marcel Hesselberth
"""


from cnav.constants import AS2RAD, SPD, MJD0, JD2000
import numpy as np
from cnav.cip import *
from cnav.dtmath import JD, MJD, TJC, TF
from cnav.webdata import leapseconds, finals


np.set_printoptions(precision = 20)


#  given:
utc_date = (2006, 1, 15)  # Wallace & Capitaine 2006
utc_time = (21, 24, 37.5)
ut1mutc = 0.3341

# calculate date / times:
TTmTAI = 32.184
TAImUTC_2006 = 33
utc_jdate = JD(*utc_date)
utc_mjd = MJD(*utc_date)
utc_time = TF(*utc_time)
tai_time = utc_time + TAImUTC_2006 / SPD
tt_time = tai_time + TTmTAI / SPD


def test_leapseconds():
    TAImUTC = leapseconds()
    ls2006 = TAImUTC(*utc_date)
    assert(ls2006 == TAImUTC_2006)
    
def test_dUT1():
    dUT1 = finals()
    dt = dUT1(utc_mjd)
    assert(abs(dt == ut1mutc) < 1e-4)  # below least significant digit

def test_date():
    ref = 53750
    assert(utc_mjd == ref)  # date exact

def test_time():
    ref = .892855138888889
    assert(abs(tt_time - ref) * SPD < 1e-10)  # time error well below 1ns

tt_mjd = utc_mjd + tt_time
tt_tjc = TJC(MJD0, tt_mjd)
def test_mjd_tjc():
    ref = 53750.892855138888889
    assert(abs(tt_mjd - ref) * SPD < 1e-10)  # time error well below 1ns
    ref = 0.06040774415164651  # Jcy
    assert(abs(tt_tjc - ref) < 1e-17)  # century

dUT1 = finals()
dt = dUT1(utc_mjd)
ut1_mjd = utc_mjd + (utc_time + ut1mutc / SPD)
def test_ut1():
    ref = 53750.892104561342593
    assert(abs(ut1_mjd - ref) * SPD < 1e-10)  # time error well below 1ns

# # continue with reference ut1mutc
ut1_2000 = (ut1mutc / SPD + utc_time, utc_jdate - JD2000)

def test_X():
    print("XYZs")
    X, Y, s = XYs(tt_tjc)
    print(f"X = {X/AS2RAD}")
    ref = +120.635997299064
    X_err = abs(X/AS2RAD - ref)
    print(f"X ERR: {X_err:.1e} arcsec\n")
    assert(X_err < 1e-7)

def test_Y():
    X, Y, s = XYs(tt_tjc)
    print(f"Y = {X/AS2RAD}")
    ref = +8.567258740044 
    Y_err = abs(Y/AS2RAD - ref)
    print(f"Y ERR: {Y_err:.1e} arcsec\n")
    assert(Y_err < 1e-7)

def test_Z():
    X, Y, s = XYs(tt_tjc)
    Z = np.sqrt(1 - X*X - Y*Y)
    print(f"Z= {Z}")
    ref = +0.999999828106893 
    Z_err = abs(Z - ref)
    print(f"Z ERR: {Z_err:.1e} (-)\n")
    assert(Z_err < 1e-7)

def test_s():
    X, Y, s = XYs(tt_tjc)
    print(f"s= {s/AS2RAD}")
    ref = -0.002571986
    s_err = abs(s/AS2RAD - ref)
    print(f"s ERR: {s_err:.1e} arcsec\n")
    assert(s_err < 1e-7)

def test_Mcio():
    print("Mcio")
    mcio = Mcio(tt_tjc)
    print(mcio)
    ref = [[+0.99999982896948063, +0.00000000032319161, -0.00058485982037403],
           [-0.00000002461548575, +0.99999999913741183, -0.00004153523474454],
           [+0.00058485981985612, +0.00004153524203735, +0.99999982810689262]]
    err = abs(mcio - ref)
    print(err, '\n')
    print("ERR Mcio")
    assert( np.max(abs(err)) < 3e-12 )
    
def test_ERA():
    print("ERA")
    #ut1_2000 = (ut1mutc / SPD + utc_time) + ()
    era = ERA(ut1_2000)
    print(f"ERA = {era/DEG2RAD} degrees")
    ref = 76.265431053522
    s_err = abs(era/DEG2RAD - ref)
    print(f"ERA ERR: {s_err:.1e} degree\n")
    assert(s_err < 1e-10)

def test_R():
    print("R")
    r = R(ut1_2000, tt_tjc)
    print(r)
    ref = [[+0.23742421473053985, +0.97140604802742432, -0.00017920749958268],
           [-0.97140588849284706, +0.23742427873021974, +0.00055827489403210],
           [+0.00058485981985612, +0.00004153524203735, +0.99999982810689262]]
    print("ERR R")
    err = abs(r - ref)
    print(r - ref, '\n')
    assert(np.max(np.abs(err)) < 3e-12)

def test_PFW():
    print("PFW")
    gamma, phi  = PFW(tt_tjc)
    print(f"gamma = {gamma/AS2RAD} arcsec")
    ref = +0.586558662
    gamma_err = abs(gamma/AS2RAD - ref)
    print(f"gamma ERR: {gamma_err:.1e} arcsec\n")
    assert(gamma_err < 1e-7)

    print(f"phi = {phi/AS2RAD} arcsec")
    ref = +84378.585257806
    phi_err = abs(phi/AS2RAD - ref)
    print(f"phi ERR: {phi_err:.1e} arcsec\n")
    assert(phi_err < 1e-7)

def test_Mclass():
    print("Mclass")
    mclass, eo = Mclass_EO(tt_tjc)
    print(mclass)
    ref = [[+0.99999892304984912, -0.00134606988972260, -0.00058480338056834],
           [+0.00134604536839225, +0.99999909318492665, -0.00004232245992880],
           [+0.00058485981924879, +0.00004153524246778, +0.99999982810689296]]
    print("ERR Mclass")
    err = mclass - ref
    np.set_printoptions(precision = 1)
    print(mclass - ref, '\n')
    np.set_printoptions(precision = 20)
    assert(np.max(np.abs(err)) < 3e-12)

def test_Mclass_Mcio():
    print("Mclass(Mcio)")
    mcio = Mcio(tt_tjc)
    mclass = Mclass(mcio, tt_tjc)
    print(mclass)
    ref = [[+0.99999892304984912, -0.00134606988972260, -0.00058480338056834],
           [+0.00134604536839225, +0.99999909318492665, -0.00004232245992880],
           [+0.00058485981924879, +0.00004153524246778, +0.99999982810689296]]
    print("ERR Mclass")
    err = mclass - ref
    print(mclass - ref, '\n')
    assert(np.max(np.abs(err)) < 3e-12)

def test_EO():
    print("EO")
    mclass, eo = Mclass_EO(tt_tjc)
    print(f"EO = {eo/AS2RAD} arcsec")
    ref = -277.646995746
    eo_err = abs(eo/AS2RAD - ref)
    print(f"EO ERR: {eo_err:.1e} arcsec\n")
    assert(eo_err < 1e-7)

def test_GST():
    print("GST")
    mclass, eo = Mclass_EO(tt_tjc)
    gst = GST(ut1_2000, eo) / DEG2RAD
    gst /= 15
    print(gst)
    hh, mm = divmod(gst, 1)
    mm, ss = divmod(mm*60, 1)
    ss*=60
    hh = int(hh)
    mm = int(mm)
    print(f"GST = {hh}:{mm}:{ss}")
    hh_ref, mm_ref, ss_ref = 5,  5, 22.213252562
    print(f"GST ERR: {hh-hh_ref}:{mm-mm_ref}:{ss-ss_ref}\n")
    assert(hh == hh_ref and mm == mm_ref and (ss - ss_ref) < 1e-7)

