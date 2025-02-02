#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 22:26:04 2025

@author: Marcel Hesselberth
"""


from cnav.dtmath import *
import pytest


def test_is_gregorian():
    assert(is_gregorian(-5000, 1, 1) is False)
    assert(is_gregorian(-4712, 1, 1) is False)
    assert(is_gregorian(0, 12, 25) is False)
    assert(is_gregorian(1000, 12, 31) is False)
    assert(is_gregorian(1582, 1, 1) is False)
    assert(is_gregorian(1582, 6, 15) is False)
    assert(is_gregorian(1582, 10, 1) is False)
    assert(is_gregorian(1582, 10, 4) is False)
    assert(is_gregorian(1582, 10, 15) is True)
    assert(is_gregorian(1582, 10, 17) is True)
    assert(is_gregorian(1582, 11, 4) is True)
    assert(is_gregorian(1582, 12, 31) is True)
    assert(is_gregorian(1602, 10, 15) is True)
    assert(is_gregorian(1941, 6, 15) is True)
    assert(is_gregorian(2025, 1, 11) is True)

def test_JD():
    assert(JD(2025,   1, 11)   == 2460686.5 )
    assert(JD(2000,   1, 1.5)  == 2451545   )
    assert(JD(1999,   1, 1)    == 2451179.5 )
    assert(JD(1987,   1, 27)   == 2446822.5 )
    assert(JD(1987,   6, 19.5) == 2446966   )
    assert(JD(1988,   1, 27)   == 2447187.5 )
    assert(JD(1988,   6, 19.5) == 2447332   )
    assert(JD(1900,   1, 1)    == 2415020.5 )
    assert(JD(1600,   1, 1)    == 2305447.5 )
    assert(JD(1600,  12, 31)   == 2305812.5 )
    assert(JD(837,    4, 10.3) == 2026871.8 )
    assert(JD(-123,  12, 31)   == 1676496.5 )
    assert(JD(-122,   1, 1)    == 1676497.5 )
    assert(JD(-1000,  7, 12.5) == 1356001.0 )
    assert(JD(-1000,  2, 29)   == 1355866.5 )
    assert(JD(-1001,  8, 17.9) == 1355671.4 )
    assert(JD(-4712,  1, 1.5)  == 0.0       )

def test_MJD():
    assert(MJD(1858, 11, 17)   == 0 )
    with pytest.raises(AssertionError) as excinfo:
        skip = JD(1582, 10, 10)

def test_JD_RJD_invariant():
    for i in range(1721057, 2634166):  # year 0 until 2500
        jd = i + 0.5
        YYYY, MM, DD, F = RJD(jd)
        jd2 = JD(YYYY, MM, DD)
        assert(jd==jd2)
        assert(F == 0)

def test_calendar_reform():
        before = JD(1582, 10, 4 )
        after  = JD(1582, 10, 15)
        for DD in range(5, 15):
            with pytest.raises(AssertionError) as excinfo:
                skip = JD(1582, 10, DD)

def test_RMJD():
    assert(RMJD(0) == (1858, 11, 17))
    
def test_TJC():
    tjc = TJC(MJD0, 53750.8928551388888888)
    ref = 0.06040774415164651
    err = abs(tjc - ref)
    assert(err < 1e-7)

def test_BY():
    jd = 50 * 365.242198781 + 2415020.31352
    by = BY(jd)
    ref = 1950
    err = abs(by - ref)
    assert(err < 1e-7)

def test_TF():
    assert ( TF(24, 0, 0) == 1 )

def test_weekday_nr():
    jd = JD(2025, 1, 12)
    assert( weekday_nr(jd) == 0 )

def test_weekday_str():
    jd = JD(2025, 1, 11)
    assert( weekday_str(jd) == "Saturday" )

def test_bisect(): 
    l = [-7, -1, 2, 4, 9, 16, 25, 32, 64, 128, 255]
    assert( bisect(l, -10) == None )
    assert( bisect(l, -7 ) == -7 )
    assert( bisect(l, 255 ) == 255 )
    assert( bisect(l, 19 ) == 16 )
    assert( bisect(l, 30 ) == 25 )
    assert( bisect(l, 45 ) == 32 )
    assert( bisect(l, 17 ) == 16 )

    
    
    