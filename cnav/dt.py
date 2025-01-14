#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:13:23 2025

@author: Marcel Hesselberth
"""

UT1 = 0
UTC = 1
TAI = 2
TT  = 3
GPS = 4
TDB = 5
TCG = 6
TCB = 7

ts = {"UT1": UT1, "UTC": UTC, "TAI": TAI, "TT":TT, "GPS":GPS, "TDB":TDB,
      "TCG":TCG, "TCB":TCB}

timescales = ts.keys()

from cnav.dtmath import MJD, RMJD

class Date(tuple):  # Date is immutable

    def __new__ (cls, yyyy, mm, dd):
        mjd = MJD(yyyy, mm, dd)
        return super(Date, cls).__new__(cls, (mjd,))

    def __init__(self, yyyy, mm, dd):
        pass



    def __str__(self):
        yyyy, mm, dd, f = RMJD(self[0])
        return f"{yyyy:= 05d}-{mm:02d}-{dd:02d}"

if __name__ == '__main__':
    date = Date(2000, 5, 1)
    print(date, date[0])
