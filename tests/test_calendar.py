#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 22:26:04 2025

@author: Marcel Hesselberth
"""


from cnav.calendar import *
import pytest


def test_is_gregorian():
    cal = Calendar()
    assert(cal.is_gregorian(-5000, 1, 1) is False)
    assert(cal.is_gregorian(-4712, 1, 1) is False)
    assert(cal.is_gregorian(0, 12, 25) is False)
    assert(cal.is_gregorian(1000, 12, 31) is False)
    assert(cal.is_gregorian(1582, 1, 1) is False)
    assert(cal.is_gregorian(1582, 6, 15) is False)
    assert(cal.is_gregorian(1582, 10, 1) is False)
    assert(cal.is_gregorian(1582, 10, 4) is False)
    assert(cal.is_gregorian(1582, 10, 15) is True)
    assert(cal.is_gregorian(1582, 10, 17) is True)
    assert(cal.is_gregorian(1582, 11, 4) is True)
    assert(cal.is_gregorian(1582, 12, 31) is True)
    assert(cal.is_gregorian(1602, 10, 15) is True)
    assert(cal.is_gregorian(1941, 6, 15) is True)
    assert(cal.is_gregorian(2025, 1, 11) is True)

def test_JD():
    cal = Calendar()
    
    assert(cal.JD(2025,   1, 11)   == 2460687, True )
    assert(cal.JD(2000,   1, 1)  == 2451545, True   )
    assert(cal.JD(1999,   1, 1)    == 2451180, True )
    assert(cal.JD(1987,   1, 27)   == 2446823, True )
    assert(cal.JD(1987,   6, 19) == 2446966, True   )
    assert(cal.JD(1988,   1, 27)   == 2447188, True )
    assert(cal.JD(1988,   6, 19) == 2447332, True   )
    assert(cal.JD(1900,   1, 1)    == 2415021, True )
    assert(cal.JD(1600,   1, 1)    == 2305448, True)
    assert(cal.JD(1600,  12, 31)   == 2305813, True)
    assert(cal.JD(837,    4, 10) == 2026872, False)
    assert(cal.JD(-123,  12, 31)   == 1676497, False)
    assert(cal.JD(-122,   1, 1)    == 1676498, False)
    assert(cal.JD(-1000,  7, 12) == 1356001, False)
    assert(cal.JD(-1000,  2, 29)   == 1355867, False)
    assert(cal.JD(-1001,  8, 18) == 1355672, False)
    assert(cal.JD(-4712,  1, 1)  == 0, False)

def test_MJD():
    cal = Calendar()
    assert(cal.MJD(1858, 11, 17)   == 0 )
    with pytest.raises(ValueError) as excinfo:
        skip = cal.JD(1582, 10, 10)

def test_JD_RJD_invariant():
    cal = Calendar()
    for i in range(1721057, 2634166):  # year 0 until 2500
        jd = i
        year, month, day = cal.RJD(jd)
        jd2 = cal.JD(year, month, day)[0]
        assert(jd==jd2)

def test_calendar_reform():
    cal = Calendar()
    before = cal.JD(1582, 10, 4 )
    after  = cal.JD(1582, 10, 15)
    for DD in range(5, 15):
        with pytest.raises(ValueError) as excinfo:
            skip = cal.JD(1582, 10, DD)

def test_RMJD():
    cal = Calendar()
    assert(cal.RMJD(0) == (1858, 11, 17))
    
# def test_TJC():
#     cal = Calendar()
#     tjc = cal.TJC(MJD0, 53750.8928551388888888)
#     ref = 0.06040774415164651
#     err = abs(tjc - ref)
#     assert(err < 1e-7)

# def test_BY():
#     cal = Calendar()
#     jd = 50 * 365.242198781 + 2415020.31352
#     by = cal.BY(jd)
#     ref = 1950
#     err = abs(by - ref)
#     assert(err < 1e-7)

def test_weekday():
    cal = Calendar()
    jd = cal.JD(2025, 1, 12)[0]
    assert( cal.weekday(jd) == 0 )

def test_weekday_str():
    cal = Calendar()
    jd = cal.JD(2025, 1, 11)[0]
    assert( cal.weekday_str(jd) == "Saturday" )

def test_bisect(): 
    l = [-7, -1, 2, 4, 9, 16, 25, 32, 64, 128, 255]
    assert( bisect(l, -10) == None )
    assert( bisect(l, -7 ) == -7 )
    assert( bisect(l, 255 ) == 255 )
    assert( bisect(l, 19 ) == 16 )
    assert( bisect(l, 30 ) == 25 )
    assert( bisect(l, 45 ) == 32 )
    assert( bisect(l, 17 ) == 16 )

    
    
    