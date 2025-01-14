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

from cnav.constants import months
from cnav.dtmath import JD, MJD, RMJD, weekday_str, weekday_nr
import re


class TimeDelta(tuple):  # TODO: add formatting code (precision)
    def __new__ (cls, days, seconds):
        return super(Date, cls).__new__(cls, (days, seconds))

    def __init__ (cls, days, seconds):
        pass

    def __str__(self):
        days, seconds = self
        return f"{days} days + {seconds} seconds"  # leave out days


class Date(tuple):  # Date is immutable
    def __new__ (cls, yyyy, mm, dd):
        mjd = MJD(yyyy, mm, dd)
        return super(Date, cls).__new__(cls, (yyyy, mm, dd))

    def __init__(self, yyyy, mm, dd):
        jd = JD(yyyy, mm, dd)
        jd_zero = JD(yyyy - 1, 12, 31)
        day_of_year = int(jd - jd_zero)
        weekday_zero = weekday_nr(jd_zero)
        first_sun = (  - weekday_zero) % 7
        first_mon = (1 - weekday_zero) % 7
        week_num_U = 1 + (day_of_year - first_sun) // 7
        week_num_W = 1 + (day_of_year - first_mon) // 7
        self.strf = {}
        self.strf["A"] = weekday_str(jd)
        self.strf["a"] = self.strf["A"][:3]
        self.strf["w"] = f"{weekday_nr(jd)}"
        self.strf["d"] = f"{dd:02d}"
        self.strf["B"] = months[mm]  # TODO: check values
        self.strf["b"] = self.strf["B"][:3]
        self.strf["m"] = f"{mm:02d}"
        self.strf["Y"] = f"{yyyy:= 05d}"
        self.strf["j"] = f"{day_of_year:03d}"
        self.strf["U"] = f"{week_num_U:02d}"
        self.strf["W"] = f"{week_num_W:02d}"

    def replace(self, x):
        code = x.group(0)[-1]
        if code in self.strf:
            return(self.strf[code])
        elif code == '%':
            return '%'
        raise(SyntaxError(f"Encountered invalid % escape ({code})"))

    def format(self, fmt):
        r = re.sub("%([a-z|A-Z|%])", self.replace, fmt)
        return(r.format())

    def __str__(self):
        yyyy, mm, dd = self
        return f"{yyyy:= 05d}-{mm:02d}-{dd:02d}"

if __name__ == '__main__':
    date = Date(2025, 1, 14)
    print(repr(date))
    print(date.format("Hello %A (%a) world! (%w)%d %B(%b) %Y-%m-%d %j %U %W%%."))
