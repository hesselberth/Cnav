#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 20:57:17 2025

@author: Marcel Hesselberth
"""

import numpy as np
from .constants import MJD0, SPD, wdays
from .cnumba import cnjit


"""
Date/time math. These algorithms use Julian day numbers to compute differences
between dates. This way a continous time scale is established starting at -4712.
A julian day starts at mean noon at the Greenwich meridian.

The Gregorian calendar reform is taken into account: the day following
October 4, 1582 (Julian calendar) is October 15, 1582.

The year before +1 is defined as the year 0 (as is usually done in astronomy).

The JD algorithm follows Meeus.

"""

@cnjit(signature_or_function='boolean(i4, i4, i4)')
def is_gregorian(YYYY:int, MM:int, DD:int) -> bool:
    """
    Check if the date is in the Gregorian calendar.
    
    Assumes that year, month and day are valid.

    Parameters
    ----------
    YYYY : int
           Year from -4712 to 9999
    MM   : int
           Month (1-12)
    DD   : int
           Day between 1 and mdays (and 29-2 in a leap year)

    Returns
    -------
    Bool
        True if the date is in the Gregorian calendar, otherwise False.

    """
    return (YYYY + MM/12) * 365.25 + DD > 578140

@cnjit(signature_or_function='f8(i4, i4, f8)')
def JD(YYYY:int, MM:int, DD:float) -> float:
    """
    Julian day number.
    
    Assumes that year, month and day are valid.
    Works for all positive JD (years from -4712).


    Parameters
    ----------
    YYYY : int
           Year.
    MM   : int
           Month.
    DD   : int
           Day.

    Returns
    -------
    float
           The corresponding Julian day number.

    """
    assert (YYYY >= -4712)
    if YYYY == 1582 and MM == 10:
        assert (DD <=4 or DD >= 15)
    if MM <= 2:
        Y  = YYYY - 1
        M  = MM + 12
    else:
        Y = YYYY
        M = MM
    A = int(Y/100)
    if (YYYY + MM/12) * 365.25 + DD < 578140:
        B = 0
    else:
        B = 2 - A + int(A/4)
    return int(365.25*(Y+4716)) + int(30.6001*(M+1)) + DD + B - 1524.5

@cnjit(signature_or_function='i4(i4, i4, i4)')
def MJD(YYYY:int, MM:int, DD:int) -> int:
    """
    Modified Julian day of a valid date.

    Parameters
    ----------
    YYYY : int
           Year
    MM   : int
           Month
    DD   : int
           Day

    Returns
    -------
    int
        Modified Julian day for the geven date

    """
    assert (YYYY >= -4712)
    if YYYY == 1582 and MM == 10:
        assert (DD <=4 or DD >= 15)
    if MM <= 2:
        Y  = YYYY - 1
        M  = MM + 12
    else:
        Y = YYYY
        M = MM
    A = int(Y/100)
    if (YYYY + MM/12) * 365.25 + DD < 578140:
        B = 0
    else:
        B = 2 - A + int(A/4)
    return int(365.25*(Y+4716)) + int(30.6001*(M+1)) + DD + B - 2401525

@cnjit(signature_or_function='Tuple((i4, i4, i4, f8))(f8)')
def RJD(jd:float) -> (int, int, int, float):
    """
    Reverse Julian Day. Compute date (YYYY, MM, DD) from jd.
    Computes the time F as a day fraction as well.

    Parameters
    ----------
    jd : float
         Julian day. jd must be positive (or zero),

    Returns
    -------
    YYYY : int
           Year.
    MM   : int
           Month.
    DD   : int
           Day.
    F    : float
           Day fraction.
    """
    jd5 = jd + 0.5
    Z   = int(jd5)
    F   = jd5 - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)  # positive
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    DD = B - D - int(30.6001 * E)
    if E < 14:
        MM = E - 1
    else:
        MM = E - 13
    assert(MM >=1 and MM <= 12)
    if MM > 2:
        YYYY = C - 4716
    else:
        YYYY = C - 4715
    return YYYY, MM, DD, F

def RMJD(mjd):
    """
    Reverse Modified Julian Day. Compute date (YYYY, MM, DD) from the mjd.
    Computes the time F as a day fraction as well.

    Parameters
    ----------
    mjd : float
         Modified Julian day.

    Returns
    -------
    YYYY : int
           Year.
    MM   : int
           Month.
    DD   : int
           Day.
    F    : float
           Day fraction.
    """
    return(RJD(mjd + MJD0))

def TJC(tt, tt2=0.0):
    return ((tt - 2451545.0) + tt2) / 36525

# Besselian year, argument is full TT julian date
def BY(tt_jd):
    return 1900.0 + (tt_jd - 2415020.31352) / 365.242198781

def TF(hh, mm, ss):
    ssum = (ss + mm * 60) + hh * 3600
    return ssum / SPD

# 0 is sunday, 1 is monday etc.
def weekday_nr(jd):    # jd is a julian day number without time
    return int((jd+1.5)) % 7

# Three character day string
def weekday_str(jd):
    return wdays[weekday_nr(jd)]

# find value in L less than or equal to jd
def bisect(l, jd): 
    if l[0] > jd:       # no smaller value exists, answer undefined
        return None
    left = 0
    right = len(l)
    while right - left > 1:
        mid = (right + left) // 2
        if l[mid] <= jd:
            left = mid
        else:
            right = mid
    return l[left]
