#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 20:57:17 2025

@author: Marcel Hesselberth
"""

import numpy as np
from math import floor
from cnav.constants import MJD0, SPD, wdays, JD2000, mdays
from cnav.cnumba import cnjit
from collections import namedtuple
from operator import index as _index

MINYEAR = -4712
MAXYEAR = 9999
MJD0 = 2400001

DateTuple = namedtuple("datetuple", ["year", "month", "day"])

"""
Calendar class for datetime.
There are 5 supported calendars:
    - Proleptic Julian: Quadrennial leap years, counts days from -4712 (=-4713BC)
    - Proleptic Gregorian: Similar to Julian plus a centennial leap year rule.
    - Mixed: A julian calendar that switches to Gregorian at a given (Gregorian)
      date. Algoritmically this date must be after March 1, 200 because otherwise
      a date would occur twice.
    - Default: The mixed calendar, where the Julian-Gregorian reform date is
      October 15, 1582.
    - The ISO calendar. This calendar is an overlay over the proleptic Gregorian
      calendar. An error will occur in the Julian period so for historical use
      of iso week numbers, set a Mixed calendar with a reform date before the
      earliest date or a Julian calendar (starting at January 1, -4712).
In all cases years from -4712 until 9999 are supported. Julian days for the
Gregorian and Julian calendars are compatible but because of the different
leap years julian days have a different date. Therefore the last days in 9999
can not be converted from Julian to Gregorian as the result will be in 10.000.
"""

class Calendar:
    default = 0
    julian = 1
    gregorian = 2
    mixed = 3
    JD0 = 1721423  # 22.5 for midnight
    GD0 = JD0 + 2
    MAXYEAR = 9999

    def __init__(self):
        self.setDefault()

    def setDefault(self):
        """
        Set the Default calendar.

        The Default calendar is a mixed calendar with a Julian to Gregorian
        calendar reform at October 15, 1582.

        Returns
        -------
        None.

        """
        self.setMixed(1582, 10, 15)

    def setJulian(self):
        """
        Set the Julian calendar from -4712-1-1 until 9999-31-1.

        Returns
        -------
        None.

        """
        self._gregorian_reform_date = DateTuple(self.MAXYEAR+1, 1, 1)
        self.reform = 1, -1  # no ValueError

    def setGregorian(self):
        """
        Set the (proleptic) Gregorian calendar from -7412-1-1 until 9999-12-31.

        Returns
        -------
        None.

        """
        self._gregorian_reform_date = DateTuple(-4712, 1, 1)
        self.reform = -1, -1

    def setMixed(self, gr_year, gr_month, gr_day):
        """
        Set a mixed Julian-Gregorian calendar.

        In some countries the calendar reform was implemented (much) after
        1582. A historically correct calendar can be set by calling Mixed
        with the Gregorian reform date.

        Parameters
        ----------
        gr_year : int
            (Gregorian) year of the calendar reform.
        gr_month : int
            (Gregorian) month of the calendar reform.
        gr_day : int
            (Gregorian) day of the calendar reform.

        Returns
        -------
        None.

        """
        year, month, day = self._check_date_fields(
            gr_year, gr_month, gr_day, self.gregorian)
        assert (year, month, day) >= (200, 3, 1)  # ensure unique dates
        self._gregorian_reform_date = DateTuple(year, month, day)
        jdl = self._GD(year, month, day)
        jdh = self._JD(year, month, day)
        self.reform = jdl, jdh

    def JD(self, year, month, day):
        """
        Julian day number of the given date according to the
        currently set calendar.

        Parameters
        ----------
        year : int

        month : int

        day : int


        Raises
        ------
        ValueError
            If the date is not valid.

        Returns
        -------
        int
            The integer julian day (= at noon).
        """
        year, month, day = self._check_date_fields(year, month, day)
        if self.is_gregorian(year, month, day):
            return self._GD(year, month, day)
        jd = self._JD(year, month, day)
        jdl, jdh = self.reform
        if jd >= jdl and jd < jdh:
            raise ValueError(f"nonexistant date ({year}-{month}-{day})")
        return jd

    def RJD(self, jd):
        """
        Reverse Julian Day.

        Given a JD, compute the date according to the currently set calendar.

        Parameters
        ----------
        jd : int
            Julian day.

        Returns
        -------
        year : int

        month : int

        day : int
        """
        if jd >= self.reform[0]:
            year, month, day = self._RGD(jd)
        else:
            year, month, day = self._RJD(jd)
        year, month, day = self._check_date_fields(year, month, day)
        return (year, month, day)
    
    def MJD(self, year, month, day):
        return self.JD(year, month, day) - MJD0

    def RMJD(self, jd):
        return self.RJD(jd + MJD0)

    # must go to date
    def TJC(self, tt, tt2=0.0):
        return ((tt - 2451545.0) + tt2) / 36525
    
    def BY(self, tt_jd):  # Besselian year, argument is full TT julian date
        return 1900.0 + (tt_jd - 2415020.31352) / 365.242198781

    def is_gregorian(self, year, month, day):
        """
        Check if the date is in the Gregorian calendar.

        This function takes into account the currently set calendar reform date.

        Parameters
        ----------
        year : int

        month : int

        day : int


        Returns
        -------
        Boolean
        True if the date is in the Gregorian calendar.
        """
        return (year, month, day) >= self._gregorian_reform_date

    def is_leapyear(self, year):
        """
        Check if the date is a leap year.

        This function takes into account the currently set calendar reform date
        and whether the Julian of Gregorian leap year rule should be applied.

        Parameters
        ----------
        year : int

        month : int

        day : int


        Returns
        -------
        Boolean
        True if the date is a leap year.
        """
        if self.is_gregorian(year, 2, 29):
            return self._is_gregorian_leapyear(year)
        return self._is_julian_leapyear(year)

    def weekday(self, jd):
        """
        The weekday number of Julian day jd.

        jd is interpreted according to the current calendar setting
        (reform date).

        Parameters
        ----------
        jd : int
            Julian day

        Returns
        -------
        int
            0: Sunday
            1: Monday
            2: Tuesday
            3: Wednesday
            4: Thursday
            5: Friday
            6: Saturday
        """
        return int((jd + 1.5)) % 7

    def weekday_str(self, jd):
        """
        The weekday (string) of Julian day jd.

        jd is interpreted according to the current calendar setting
        (reform date).


        Parameters
        ----------
        jd : int
            Julian day

        Returns
        -------
        string            
            "Sunday"
            "Monday"
            "Tuesday"
            "Wednesday"
            "Thursday"
            "Friday"
            "Saturday"
        """
        return wdays[self.weekday(jd)]

    def isoweekday(self, jd):
        """
        The ISO weekday number of Julian day jd.

        jd is interpreted according to the current calendar setting
        (reform date).

        Parameters
        ----------
        jd : int
            Julian day

        Returns
        -------
        int
            1: Monday
            2: Tuesday
            3: Wednesday
            4: Thursday
            5: Friday
            6: Saturday
            7: Sunday
        """
        return int((jd + 0.5)) % 7 + 1

    def I2G(self, year, week, day):
        """
        Convert an ISO date to a (proleptic) Gregorian date.

        Parameters
        ----------
        year : int
            ISO year.
        week : int
            ISO week number.
        day : int
            ISO day.

        Raises
        ------
        TypeError
            Argument does not have int type.
        ValueError
            Not a valid ISO date, or outside the Gregoriann calendar
            (according to the current calendar settings).

        Returns
        -------
        year, month, day
            Gregorian date.
        """
        for var in year, week, day:
            if not isinstance(var, int):
                raise TypeError
        if not MINYEAR <= year <= MAXYEAR:
            raise ValueError
        if not 1 <= day <= 7:
            raise ValueError
        if not 1 <= week <= self._iso_weeks_in_year(year):
            raise ValueError
        # jd of 1-1 of the year containing thursday of week
        jd_1_1 = self._GD(year, 1, 1)
        # ISO weekday of jan_1 of year containing thursday of week
        weekday_1_1 = self.isoweekday(jd_1_1)
        # day in january of first thursday of year of thursday of week
        first_thursday = int((4 - weekday_1_1)) % 7 + 1
        # jd of first thursday of year containing thursday of week
        jd_first_thursday = self._GD(year, 1, first_thursday)
        # number of days to date, counted from 1st thursday of year
        days_from_first_thursday = 7 * (week - 1) + day - 4
        # jd of gregorian date
        jd = jd_first_thursday + days_from_first_thursday
        return self._RGD(jd)

    def G2I(self, gr_year, gr_month, gr_day):  # TODO: check values
        """
        Convert a (proleptic) Gregorian date to an ISO date.

        Parameters
        ----------
        year : int

        month : int

        day : int


        Raises
        ------
        ValueError
            Not a valid Gregorian date, or outside the Gregoriann calendar
            (according to the current calendar settings).

        Returns
        -------
        ISO year, ISO weeknum, ISO day
            ISO date.
        """
        if not self.is_gregorian(gr_year, gr_month, gr_day):
            raise ValueError("iso date not Gregorian, change calendar setting")
        year, month, date = self._check_date_fields(gr_year, gr_month, gr_day)
        jd = self._GD(gr_year, gr_month, gr_day)
        iso_weekday = self.isoweekday(jd)
        jd_thursday = jd + 4 - iso_weekday
        iso_year, m, d = self._RGD(jd_thursday)
        jd_1_1 = self._GD(iso_year, 1, 1)
        weekday_1_1 = self.isoweekday(jd_1_1)
        first_thursday = (4 - weekday_1_1) % 7 + 1
        jd_first_thursday = self._GD(iso_year, 1, first_thursday)
        iso_weeknum = (jd_thursday - jd_first_thursday) // 7 + 1
        return (iso_year, iso_weeknum, iso_weekday)

    def _check_date_fields(self, year, month, day, calendar=mixed):
        year = _index(year)
        month = _index(month)
        day = _index(day)
        if not MINYEAR <= year <= MAXYEAR:
            raise ValueError('year must be in %d..%d' %
                             (MINYEAR, MAXYEAR), year)
        if not 1 <= month <= 12:
            raise ValueError('month must be in 1..12', month)
        dmax = mdays[month]
        if calendar in [self.default, self.mixed]:
            is_leapyear = self.is_leapyear
        elif calendar == self.gregorian:
            is_leapyear = self._is_gregorian_leapyear
        elif calendar == self.julian:
            is_leapyear = self._is_julian_leapyear
        else:
            raise ValueError(
                "calendar not in [default, julian, gregorian, mixed]")
        if month == 2 and is_leapyear(year):
            dmax += 1
        if not 1 <= day <= dmax:
            raise ValueError('day must be in 1..%d' % dmax, day)
        return year, month, day

    def _iso_weeks_in_year(self, year):
        py = year + floor(year/4) - floor(year/100) + floor(year/400)
        y = year - 1
        pym1 = y + floor(y/4) - floor(y/100) + floor(y/400)
        if (py % 7 == 4) or (pym1 % 7 == 3):
            return 53
        return 52

    def _is_gregorian_leapyear(self, year):
        """
        Check if a given year is a leap year in the (proleptic)
        Gregorian calendar.

        Parameters
        ----------
        year : int
            DESCRIPTION.

        Returns
        -------
        Boolean
            True if year is a leap year
        """
        if year % 4 == 0:               # possibly leap
            if year % 400 == 0:         # leap
                leapyear = True
            else:
                if year % 100 == 0:     # common
                    leapyear = False
                else:                   # leap
                    leapyear = True
        else:
            leapyear = False
        return leapyear

    def _GD0(self, gr_year, gr_month, gr_day):
        """
        Given a Gregorian date, compute the day number.
        Handles positive and negative years.
        Day 0 is at 1-1-0.

        Parameters
        ----------
        gryear : int

        grmonth : int

        grday : int


        Returns
        -------
        int
            The day number corresponding to the Gregorian date
        """
        y = gr_year - 1
        ord = 365 * y + y // 4 - y // 100 + y // 400
        ord += (367 * gr_month - 362) // 12
        if gr_month > 2:
            if self._is_gregorian_leapyear(gr_year):
                ord -= 1
            else:
                ord -= 2
        ord += gr_day
        return ord

    def _GD(self, gr_year, gr_month, gr_day):
        return self._GD0(gr_year, gr_month, gr_day) + self.GD0

    def _RGD(self, jd):
        """
        Given a JD number, compute the corresponding Gregorian Date.

        RGD(GD(y, m, d)) is an invariant. This function can be used to
        compute a proleptic Gregorian date corresponding to a Julian date.

        Parameters
        ----------
        jd : int
             Day number.

        Returns
        -------
        Gregorian date.
        """
        d = jd - self.GD0 - 1
        n400 = d // 146097
        d1 = d % 146097
        n100 = d1 // 36524
        d2 = d1 % 36524
        n4 = d2 // 1461
        d3 = d2 % 1461
        n1 = d3 // 365
        y = 400 * n400 + 100 * n100 + 4 * n4 + n1
        if not (n1 == 4 or n100 == 4):
            y += 1
        y = int(y)
        if n1 != 4 and n100 != 4:
            p = d3 % 365
        else:
            return y, 12, 31
        if self._is_gregorian_leapyear(y):
            if p < 60:
                c = 0
            else:
                c = 1
        else:
            if p < 59:
                c = 0
            else:
                c = 2
        m = (12 * (p + c) + 373) // 367
        d = d - self._GD0(y, m, 0) + 1
        return (y, int(m), int(d))

    def date_from_gregorian(self, gr_year, gr_month, gr_day):
        return self.RJD(self.JD(self, gr_year, gr_month, gr_day))[:3]

    def date_to_gregorian(self, year, month, day):
        jd = self.JD(year, month, day)
        rgd = self.RJD(jd)
        return rgd

    # algorithms for the proleptic Julian calendar until -4712

    def _is_julian_leapyear(self, year):
        """
        Check if a given year is a leap year in the (proleptic) Julian calendar

        This funtion assumes counting of BCE years starts at zero.
        1 BCE (0) is a leap year.

        Parameters
        ----------
        year : int

        Returns
        -------
        Boolean
            True if year is a leap year
        """
        return year % 4 == 0

    def _JD0(self, j_year, j_month, j_day):
        """
        Given a Julian date, compute the day number.
        Handles negative dates.
        0 is at 1-1-0.

        Parameters
        ----------
        j_year : int

        j_grmonth : int

        j_grday : int


        Returns
        -------
        int
            The Julian day number in the Julian calendar corresponding to the
            Gregorian date
        """
        y = j_year - 1
        ord = 365 * y + y // 4
        ord += (367 * j_month - 362) // 12
        if j_month > 2:
            if self._is_julian_leapyear(j_year):
                ord -= 1
            else:
                ord -= 2
        ord += j_day
        return ord

    def _JD(self, j_year, j_month, j_day):
        return self._JD0(j_year, j_month, j_day) + self.JD0

    def _RJD(self, jd):
        """
        Given a JD number, compute the corresponding Julian Date.

        RJD(JD(y, m, d)) is an invariant. This function can be used to
        compute a proleptic Gregorian date corresponding to a Julian date.

        Parameters
        ----------
        jd : int
             Day number.

        Returns
        -------
        TYPE: 3-tuple if type int (year, month, day)
        Julian date.
        """
        d = jd - self.JD0 - 1
        year = (4 * d + 1464) // 1461
        year = int(year)
        p = d - self._JD0(year, 1, 0)
        if self._is_julian_leapyear(year):
            if p < 60:
                c = 0
            else:
                c = 1
        else:
            if p < 59:
                c = 0
            else:
                c = 2
        m = (12 * (p + c) + 373) // 367
        d = d - self._JD0(year, m, 0) + 1
        return (year, int(m), int(d))


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
