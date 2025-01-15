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

from cnav.constants import months, MJD0, SPD
from cnav.dtmath import JD, RJD, MJD, RMJD, weekday_str, weekday_nr
from functools import total_ordering
import re


class TimeDelta(tuple):  # TODO: add formatting code (precision)
    def __new__ (cls, days, fraction):
        if not (isinstance(days, int) and isinstance(fraction, (int, float))):
            raise (TypeError("TimeDelta(int, float) expected"))
        return super(TimeDelta, cls).__new__(cls, (days, fraction))

    def days(self):
        return self[1] + self[0]
    
    def seconds(self):
        return self[1] * SPD + self[0] * SPD

    def __init__ (cls, days, seconds):
        pass

    def __str__(self):
        return f"TimeDelta{tuple.__str__(self)}"

    def __add__(self, b):
        if isinstance(b, int):
            return TimeDelta(self[0]+b, self[1])
        elif isinstance(b, float):
            if abs(b) < 1:
                return(TimeDelta(self[0], self[1] + b))
            d1, f = divmod(b, 1)
            d2, f = divmod(self[1] + f, 1)
            return TimeDelta(self[0] + d1 + d2, f)
        elif isinstance(b, TimeDelta):
            return(TimeDelta(self[0] + b[0], self[1] + b[1]))
        elif isinstance(b, Date):
            if b[1] != 0:
                raise (ValueError("{b} has nonzero time in date arithmetic"))            
            return Date(*RMJD(self.mjd + b[0])[:3])
        raise (TypeError("DateTime instance can't be added to {type(b)}"))

    def __sub__(self, b):
        if isinstance(b, int):
            return TimeDelta(self[0]+b, self[1])
        elif isinstance(b, float):
            if abs(b) < 1:
                return(TimeDelta(self[0], self[1] + b))
            d1, f = divmod(b, 1)
            d2, f = divmod(self[1] + f, 1)
            return TimeDelta(self[0] + d1 + d2, f)
        elif isinstance(b, TimeDelta):
            return(TimeDelta(self[0] + b[0], self[1] + b[1]))
        elif isinstance(b, Date):
            if b[1] != 0:
                raise (ValueError("{b} has nonzero time in date arithmetic"))            
            return Date(*RMJD(self.mjd + b[0])[:3])
        raise (TypeError("DateTime instance can't be added to {type(b)}"))

    def __eq__(self, x):
        if isinstance(x, TimeDelta):
            pass

class Date(tuple):  # Date is immutable
    def __new__ (cls, yyyy, mm, dd):
        return super(Date, cls).__new__(cls, (yyyy, mm, dd))

    def __init__(self, yyyy, mm, dd):
        self.mjd = MJD(yyyy, mm, dd)
        mjd_zero = MJD(yyyy - 1, 12, 31)
        day_of_year = int(self.mjd - mjd_zero)
        weekday_zero = weekday_nr(mjd_zero + MJD0)
        first_sun = (  - weekday_zero) % 7
        first_mon = (1 - weekday_zero) % 7
        week_num_U = 1 + (day_of_year - first_sun) // 7
        week_num_W = 1 + (day_of_year - first_mon) // 7
        self.strf = {}
        self.strf["A"] = weekday_str(self.mjd + MJD0)
        self.strf["a"] = self.strf["A"][:3]
        self.strf["w"] = f"{weekday_nr(self.mjd + MJD0)}"
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
    
    def __sub__(self, b):
        if isinstance(b, int):
            return Date(*RMJD(self.mjd - b)[:3])
        elif isinstance(b, Date):
            return TimeDelta(self.mjd - b.mjd, 0)
        elif isinstance(b, TimeDelta):
            if b[1] != 0:
                raise (ValueError("{b} has nonzero time in date arithmetic"))
            return Date(*RMJD(self.mjd - b[0])[:3])

    def __add__(self, b):
        if isinstance(b, int):
            return Date(*RMJD(self.mjd + b)[:3])
        elif isinstance(b, Date):
            return TimeDelta(self.mjd + b.mjd, 0)
        elif isinstance(b, TimeDelta):
            if b[1] != 0:
                raise (ValueError("{b} has nonzero time in date arithmetic"))            
            return Date(*RMJD(self.mjd + b[0])[:3])

if __name__ == '__main__':
    date1 = Date(2025, 1, 15)
    date2 = Date(2025, 1, 17)
    print(repr(date1))
    delta = date1 - date2
    print(delta)
    print(date2 + delta, date1 + delta)
    print(delta+delta)
    #print(date1.format("Hello %A (%a) world! (%w)%d %B(%b) %Y-%m-%d %j %U %W%%."))
