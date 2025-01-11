#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 20:57:17 2025

@author: Marcel Hesselberth
"""

import numpy as np
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

def JD(YYYY:int, MM:int, DD:int) -> float:
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
           The corresponding Julian day number

    """
    if YYYY < -4712:
        raise(ValueError("Invalid date: %d-%d-%d" % (YYYY, MM, DD)))
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