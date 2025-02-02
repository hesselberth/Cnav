#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 17:37:29 2025

@author: hessel
"""

from cnav.dtmath import MJD, GD, JD, RJD, RJDg, JDg, JDj, RJDj, date_from_gregorian, date_to_gregorian, date_from_julian, date_from_julian, is_leapyear, is_gregorian_leapyear
from cnav.constants import mdays

date = (0, 12, 31)

# for year in range(-5000, 5000 + 1):
#     for month in [1, 2, 3, 5, 6, 11, 12]:
#         dmax = mdays[month]
#         if month == 2 and is_gregorian_leapyear(year):
#             dmax += 1
#         for day in range(1, dmax+1):
#             date = (year, month, day)
#             y, m, d = RJDg(JDg(*date))
#             if not (y == year and m == month and d == day):
#                 print("RJDg(JDg)", date , (y, m, d))
#             y, m, d = RJDj(JDj(*date))
#             if not (y == year and m == month and d == day):
#                 print("RJDj(JDj)", date , (y, m, d))
#     jdl = JDg(-500, 1, 1)
#     jdh = JDg(500, 1, 1)
#     jd = jdl
# while jd <= jdh:
#     calc = JDg(*RJDg(jd))
#     if not (jd == calc):
#         print("JDg(RJDg)", jd, calc, RJDg(jd), RJDg(calc))
#     calc = JDj(*RJDj(jd))
#     if not (jd == calc):
#         print("JDj(RJDj)", jd, calc, RJDj(jd), RJDj(calc))
#     jd += 1
    
print(RJDj(JDj(2000, 2, 2)))
print(date_from_julian(2025, 2, 2))