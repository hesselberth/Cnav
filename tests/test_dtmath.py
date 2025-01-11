#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 22:26:04 2025

@author: Marcel Hesselberth
"""


from cnav.dtmath import is_gregorian, JD


def test_is_gregorian():
    assert(is_gregorian(-5000, 1, 1) is False)
    assert(is_gregorian(-4712, 1, 1) is False)
    assert(is_gregorian(0, 12, 26) is False)
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
