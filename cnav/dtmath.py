#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 20:57:17 2025

@author: hessel
"""

import numpy as np


"""
Date/time math. These algorithms use Julian day numbers to compute differences
between dates. This way a continous time scale is established starting at -4712.
A julian day starts at noon at the Greenwich meridian.

The Gregorian calendar reform is taken into account between October 4 and
October 15 1582.

The year before +1 is defined as the year 0 (as is usually done in astronomy).

The JD algorithm follows Meeus.

"""


