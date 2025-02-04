#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:13:23 2025

@author: Marcel Hesselberth
"""

from cnav.constants import months, MJD0, SPD, mdays
from cnav.dtmath import JD, RJD, MJD, RMJD, weekday_str, weekday_nr, is_gregorian, date_from_gregorian, Calendar
import re
import numpy as np
from functools import total_ordering, cached_property
import time
from collections import namedtuple
import datetime
from locale import nl_langinfo, D_T_FMT, D_FMT, T_FMT
from decimal import Decimal, InvalidOperation
from enum import Enum
from datetime import timezone, timedelta
from operator import index as _index
from math import floor
import _strptime

class ts(Enum):
    UT1 = 0
    UTC = 1
    TAI = 2
    TT = 3
    GPS = 4
    TDB = 5
    TCG = 6
    TCB = 7


tsname = {"UT1": ts.UT1, "UTC": ts.UTC, "TAI": ts.TAI, "TT": ts.TT,
          "GPS": ts.GPS, "TDB": ts.TDB, "TCG": ts.TCG, "TCB": ts.TCB}

timescales = tuple(tsname.keys())

MJD_UNIX_EPOCH = 40587  # 1-1-1970
DELTA_MAX_DAYS = 5373485

datetime.MINYEAR = -4712
MINYEAR = -4712
MAXYEAR = 9999

# IsoCalendarDate = namedtuple("IsoCalendarDate", ["year", "week", "weekday"])
struct_time = namedtuple("struct_time", ["tm_year", "tm_mon", "tm_mday",
                                         "tm_hour", "tm_min", "tm_sec",
                                         "tm_wday", "tm_yday", "tm_isdst"])

class IsoCalendarDate(tuple):
    def __new__(cls, year, week, weekday, /):
        return super().__new__(cls, (year, week, weekday))

    @property
    def year(self):
        return self[0]

    @property
    def week(self):
        return self[1]

    @property
    def weekday(self):
        return self[2]

    def __reduce__(self):
        return (tuple, (tuple(self),))

    def __repr__(self):
        return (f'{self.__class__.__name__}'
                f'(year={self[0]}, week={self[1]}, weekday={self[2]})')

class TZInfo:
    pass


_USEC = 1000  # ns
_MSEC = 1000 * _USEC
_SEC = 1000 * _MSEC
_MIN = 60 * _SEC
_HR = 60 * _MIN
_DAY = 24 * _HR
_WEEK = 7 * _DAY


class DTMeta(type):
    def __init__(cls, name, bases, dct):
        cls.cname = name
        cls.mname = cls.__module__
        super().__init__(name, bases, dct)

def _get_class_module(self):
    return self.__class__.__module__

def _check_date_fields(year, month, day):
    year = _index(year)
    month = _index(month)
    day = _index(day)
    if not MINYEAR <= year <= MAXYEAR:
        raise ValueError('year must be in %d..%d' % (MINYEAR, MAXYEAR), year)
    if not 1 <= month <= 12:
        raise ValueError('month must be in 1..12', month)
    dim = mdays[month]
    if month == 2 and is_leapyear(year):
        dim += 1
    if not 1 <= day <= dim:
        raise ValueError('day must be in 1..%d' % dim, day)
    if year == 1582 and month == 10 and day > 4 and day < 15:
        raise (ValueError("Invalid date (Gregorian reform)"))
    return year, month, day

def _check_time_fields(hour, minute, second, nanosecond, fold):
    hour = _index(hour)
    minute = _index(minute)
    second = _index(second)
    nanosecond = _index(nanosecond)
    if not 0 <= hour <= 23:
        raise ValueError('hour must be in 0..23', hour)
    if not 0 <= minute <= 59:
        raise ValueError('minute must be in 0..59', minute)
    if hour == 23 and minute == 59:
        maxsec = 60
    else:
        maxsec = 59
    if not 0 <= second <= maxsec:
        raise ValueError(f'second must be in 0..{maxsec}', second)
    if not 0 <= nanosecond <= 999999999:
        raise ValueError('microsecond must be in 0..999999999', nanosecond)
    if fold not in (0, 1):
        raise ValueError('fold must be either 0 or 1', fold)
    return hour, minute, second, nanosecond, fold

def _check_tzinfo_arg(tz):
    if tz is not None and not isinstance(tz, TZInfo):
        raise TypeError("tzinfo argument must be None or of a tzinfo subclass")
def _p(y): return y + floor(y/4) - floor(y/100) + floor(y/400)

def iso_weeks_in_year(year):
    if ( _p(year) % 7 == 4 ) or ( _p(year-1) % 7 == 3 ):
        return 53
    return 52
    
@total_ordering  # TODO: check efficiency
class TimeDelta(datetime.timedelta, metaclass=DTMeta):
    #__slots__ = "_days", "_seconds", "_microseconds", "_total_nanoseconds", "_hashcode"

    def __new__(cls, days=0, seconds=0, microseconds=0, nanoseconds=0,
                milliseconds=0, minutes=0, hours=0, weeks=0):
        
        for name, value in (
            ("days", days),
            ("seconds", seconds),
            ("microseconds", microseconds),
            ("milliseconds", milliseconds),
            ("minutes", minutes),
            ("hours", hours),
            ("weeks", weeks),
            ("nanoseconds", nanoseconds)
            ):
            if not isinstance(value, (int, float, Decimal)):
                raise TypeError(
                    f"unsupported type for timedelta {name} argument: {type(value).__name__}"
                )

        ns = Decimal()
        if nanoseconds != 0:
            ns += Decimal(nanoseconds)
        if microseconds != 0:
            ns += Decimal(microseconds) * _USEC
        if milliseconds != 0:
            ns += Decimal(milliseconds) * _MSEC
        if seconds != 0:
            ns += Decimal(seconds) * _SEC
        if minutes != 0:
            ns += Decimal(minutes) * _MIN
        if hours != 0:
            ns += Decimal(hours) * _HR
        if days != 0:
            ns += Decimal(days) * _DAY
        if weeks != 0:
            ns += Decimal(weeks) * _WEEK

        total_nanoseconds = round(ns)
        days, ns = divmod(total_nanoseconds, _DAY)
        seconds, ns = divmod(ns, _SEC)
        microseconds, ns = divmod(ns, _USEC)
        
        if abs(days) > 999999999:
            raise OverflowError(f"timedelta # of days is too large: {days}")

        self = datetime.timedelta.__new__(cls, days=days, seconds=seconds, microseconds=microseconds)
        self._total_nanoseconds = total_nanoseconds
        self._hashcode = -1
        return self

    def total_days(self):
        return Decimal.__truediv__(self, _DAY)

    def __str__(self):
        l = []
        ns =  int(self._total_nanoseconds % 1000000000)
        hours, seconds = divmod(self.seconds, 3600)
        minutes, seconds = divmod(seconds, 60)
        if self.days:
            if self.days in [-1, 1]:
                l.append(f"{self.days} day, ")
            else:
                l.append(f"{self.days} days, ")
        l.append(f"{hours}:{minutes:02d}:{seconds:02d}")
        if ns:
            ns = f".{ns:09d}"
            l.append(ns)
        return "".join(l)

    def __repr__(self):
        f = []
        #weeks, days = divmod(self.days, 7)
        days = self.days
        seconds = self.seconds
        microseconds = self.microseconds
        nanoseconds = self._total_nanoseconds % 10**3
        if days:
            f.append(f"days={days}")
        if seconds:
            f.append(f"seconds={seconds}")
        if microseconds:
            f.append(f"microseconds={microseconds}")
        if nanoseconds:
            f.append(f"nanoseconds={nanoseconds}")
        if not f:
            f.append("0")
        return f"{self.mname}.{self.cname}({', '.join(f)})"

    @classmethod
    def total_nanoseconds(self, td=None):
        if td is None:
            return self._total_nanoseconds
        if isinstance(td, TimeDelta):
            return td._total_nanoseconds
        if isinstance(td, datetime.timedelta):
            us = Decimal(td.microseconds)
            us += td.seconds * 10**6
            us += td.days * SPD * 10**6
            return 1000 * us
        raise NotImplementedError
    
    def __add__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            ns = self._total_nanoseconds + self.total_nanoseconds(other)
            return TimeDelta(nanoseconds=ns)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            ns = self._total_nanoseconds - self.total_nanoseconds(other)
            return TimeDelta(nanoseconds=ns)
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            return -self + other
        return NotImplemented


    def __neg__(self):
        return TimeDelta(nanoseconds=-self._total_nanoseconds)

    # pos inherited
    # abs inherited

    def __mul__(self, other):
        if isinstance(other, (int, float, Decimal)):
            ns = self._total_nanoseconds * other
            return TimeDelta(nanoseconds=ns)
        return NotImplemented

    __rmul__ = __mul__

    def __floordiv__(self, other):
        if isinstance(other, int):
            return TimeDelta(nanoseconds=self._total_nanoseconds // other)
        if isinstance(other, (TimeDelta, timedelta)):
            return int(self._total_nanoseconds // self.total_nanoseconds(other))
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (int, float, Decimal)):
            ns = self._total_nanoseconds / other
            return TimeDelta(nanoseconds=ns)
        if isinstance(other, (TimeDelta, timedelta)):
            return self._total_nanoseconds / self.total_nanoseconds(other)
        return NotImplemented

    # Beware: decimal can have negative remainders after mod!
    def __mod__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            return TimeDelta(nanoseconds = int(self._total_nanoseconds) % int(self.total_nanoseconds(other)))
        raise (TypeError(f"Incompatible types for % : <{self.cname}>, {type(other)}"))

    def __divmod__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            div, mod = divmod(int(self._total_nanoseconds), int(self.total_nanoseconds(other)))
            return int(div), TimeDelta(nanoseconds=mod)
        raise (TypeError(f"Incompatible types for divmod : <{self.cname}>, {type(other)}"))

    def __float__(self):
        return float(Decimal.__truediv__(self._total_nanoseconds, _DAY))

    def __eq__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            return self._total_nanoseconds == self.total_nanoseconds(other)
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (TimeDelta, timedelta)):
            return self._total_nanoseconds > self.total_nanoseconds(other)
        return NotImplemented

    def __bool__(self):
            return self._total_nanoseconds != 0

    def _getstate(self):
        return (0, 0, 0, self._total_nanoseconds)

    def __reduce__(self):
        return (self.__class__, self._getstate())
    
    def __hash__(self):
        if self._hashcode == -1:
            self._hashcode = hash(self._getstate())
        return self._hashcode

TimeDelta.min = TimeDelta(-999999999)
TimeDelta.max = TimeDelta(days=999999999, hours=23, minutes=59, seconds=59,
                          nanoseconds=999999999)
TimeDelta.resolution = TimeDelta(nanoseconds=1)

class DateTime:
    pass

@total_ordering
class Date(metaclass = DTMeta):
    #__slots__ = '_year', '_month', '_day', '_hashcode', "cname", "mjd"
    
    resolution = TimeDelta(days=1)
    ORD0 = MJD(1, 1, 2)
    dpat = "([+-]?)([0-9]{4})(?P<dash>-?)(W?)([0-5][0-9])(?:(?P=dash)([0-9]{1,2}))?"
    datepattern = re.compile(dpat)
    min_mjd = MJD(MINYEAR, 1, 1)
    max_mjd = MJD(MAXYEAR, 12, 31)
    valid_mjd = range(min_mjd, max_mjd+1)

    # def __new__(cls, year=2000, month=1, day=1):
    #     return super().__new__(cls)

    def __new__(cls, year, month=None, day=None, *args):  # args for str
        if month is None and day is None:
            if isinstance(year, int):
                self = object.__new__(cls)
                self.__setstate(year)
                self._hashcode = -1
                return self
        year, month, day = _check_date_fields(year, month, day)
        self = object.__new__(cls)
        self.year = year
        self.month = month
        self.day = day
        self._hashcode = -1
        return self

    @cached_property
    def mjd(self):
        return MJD(self.year, self.month, self.day)

    @classmethod
    def today(self):
        return self.fromtimestamp(time.time())

    @classmethod
    def fromtimestamp(self, timestamp):
        if timestamp is None:
            raise TypeError("'NoneType' object cannot be interpreted as a time stamp")
        y, m, d, hh, mm, ss, weekday, jday, dst = time.localtime(timestamp)
        return Date(y, m, d)

    @classmethod
    def fromjd(self, jd):
        return Date(*RJD(jd)[:3])

    @classmethod
    def frommjd(cls, mjd):
        if not mjd in cls.valid_mjd:
            raise ValueError("Date ordinal out of range")
        return cls(*RMJD(mjd))

    @classmethod
    def fromordinal(cls, ordinal):
        return cls.frommjd(ordinal)

    @classmethod
    def fromGregorian(cls, year, month, day):
        date = date_from_gregorian(year, month, day)
        return cls(*date)

    def jd(self):
        print(self, self.year, self.month, self.day)
        return JD(self.year, self.month, self.day)

    @classmethod
    def fromisoformat(cls, date_string):
        m = cls.datepattern.match(date_string)
        if m:
            sign, year, dash, weekformat, wm, day = m.groups()
            year = int(year) * (1 - 2 * (year == '-'))
            if year < MINYEAR or year > MAXYEAR:
                raise(ValueError("Year out of range [-4712, +9999]"))
            if weekformat == 'W':
                if day == None:
                    day = "1"
                week, day  = int(wm), int(day)
                last = Date(year, 12, 28).isocalendar().week
                if week < 1 or week > last or day < 1 or day > 7:
                    raise(ValueError(f"Date not valid ({date_string})"))
                return cls.fromisocalendar(year, week, day)
            if day == None:
                raise ValueError(f"Date not valid (day missing in '{date_string}')")
            month, day = int(wm), int(day)
            return cls(year, month, day)
        raise (ValueError(f"Error parsing ISO date {date_string}"))

    @classmethod
    def fromisocalendar(self, year, week, day=1):
        for var in year, week, day:
            if not isinstance(var, int):
                raise TypeError
        if (year < 1582) or ((year == 1582) and (week < 40 or (week == 40 and day < 5))):
            raise ValueError("ISO calendar date must be Gregorian (>= 1582-W40-5)") 
        if year > 9999:
            raise ValueError
        if week < 1 or week > iso_weeks_in_year(year) or day < 1 or day > 7:
            raise ValueError            
        thursday_year = year
        weekday_1 = Date(year, 1, 1).isoweekday()
        first_thursday = (4 - weekday_1) % 7 + 1
        first_thursday_date = Date(year, 1, first_thursday)
        return first_thursday_date + TimeDelta(7 * (week - 1) + day - 4)

    @classmethod
    def strptime(cls, date_string, format):
        args = _strptime._strptime(date_string, format)[0][:3]
        return cls(*args)

    @cached_property
    def strf(self):
        mjd_zero = MJD(self.year - 1, 12, 31)
        self.yday = int(self.mjd - mjd_zero)
        weekday_zero = weekday_nr(mjd_zero + MJD0)
        first_sun = (- weekday_zero) % 7
        first_mon = (1 - weekday_zero) % 7
        week_num_U = 1 + (self.yday - first_sun) // 7
        week_num_W = 1 + (self.yday - first_mon) // 7
        strf = {}
        strf["A"] = weekday_str(self.mjd + MJD0)
        strf["a"] = strf["A"][:3]
        strf["w"] = f"{weekday_nr(self.mjd + MJD0)}"
        strf["d"] = f"{self.day:02d}"
        strf["e"] = f"{self.day: 2d}"
        strf["B"] = months[self.month]  # TODO: check values
        strf["b"] = strf["B"][:3]
        strf["m"] = f"{self.month:02d}"
        strf["Y"] = f"{self.year:=-04d}"
        strf["y"] = strf["Y"][-2:]
        strf["H"] = "00"
        strf["I"] = "00"
        strf["p"] = "AM"
        strf["M"] = "00"
        strf["S"] = "00"
        strf["f"] = "0" * 9
        strf["z"] = ""
        strf["j"] = f"{self.yday:03d}"
        strf["U"] = f"{week_num_U:02d}"
        strf["W"] = f"{week_num_W:02d}"
        strf["z"] = ""
        strf["Z"] = ""
        strf["F"] = f"{self.year:=-04d}-{self.month:02d}-{self.day:02d}"
        strf["C"] = f"{self.year // 100:02d}"
        return strf

    def repl(self, escape):  # replace format % escape by value
        code = escape.group(0)[-1]
        if code in self.strf:
            return (self.strf[code])
        elif code == '%':
            return '%'

        # iso, not in dict because these trigger an exception for dates in
        # the Julian calendar
        elif code == "G":
            return f"{self._iso[0]:-04d}"
        elif code == "u":
            return f"{self._iso[1]:1d}"
        elif code == "V":
            return f"{self._iso[2]:02d}"

        # there is a theoretical possibility that locale could cause
        # infinite recursion here.
        elif code == 'c':
            return self.strftime(nl_langinfo(D_T_FMT))
        elif code == 'x':
            return self.strftime(nl_langinfo(D_FMT))
        elif code == 'X':
            return self.strftime(nl_langinfo(T_FMT))
        raise (SyntaxError(f"Encountered invalid % escape ({code})"))

    def strftime(self, format=None):
        if format is None:
            raise TypeError("Call to strftime without argument")
        r = re.sub("%:?([a-z|A-Z|%])", self.repl, format)
        return (r.format())

    def replace(self, **kwargs):
        year = kwargs["year"] if "year" in kwargs else self.year
        month = kwargs["month"] if "month" in kwargs else self.month
        day = kwargs["day"] if "day" in kwargs else self.day
        return type(self)(year, month, day)

    __replace__ = replace

    def timetuple(self):
        mjd_zero = MJD(self.year - 1, 12, 31)
        self.yday = self.mjd - MJD(self.year - 1, 12, 31)
        return struct_time(self.year, self.month, self.day, 0, 0, 0,
                           self.weekday(), self.yday, -1)

    def toordinal(self):
        return self.mjd

    def weekday(self):
        return (weekday_nr(self.mjd + MJD0) - 1) % 7

    def isoweekday(self):
        return (weekday_nr(self.mjd + MJD0) - 1) % 7 + 1

    def isgregorian(self):
        return is_gregorian(self.year, self.month, self.day)

    @cached_property
    def _iso(self):
        if not self.isgregorian():
            raise ValueError("ISO notation is only defined for Gregorian dates (>= 1582-10-15)")
        iso_weekday = (weekday_nr(self.mjd + MJD0) - 1) % 7 + 1
        thursday_date = self + TimeDelta(4 - iso_weekday)
        iso_year = thursday_date.year
        weekday_1 = Date(iso_year, 1, 1).isoweekday()
        first_thursday = (4 - weekday_1) % 7 + 1
        first_thursday_date = Date(iso_year, 1, first_thursday)
        iso_weeknum = (thursday_date - first_thursday_date).days // 7 + 1
        return IsoCalendarDate(iso_year, iso_weeknum, iso_weekday)

    def isocalendar(self):
        return self._iso

    def isoformat(self):
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"

    def ctime(self):
        return f"{self.strf['a']} {self.strf['b']} {self.day:2d} 00:00:00 {self.year}"

    def __str__(self):
        if is_gregorian(self.year, self.month, self.day):
            cal = "G"
        else:
            cal = "J"
        return f"{cal}{self.year: 04d}-{self.month:02d}-{self.day:02d}"

    def __repr__(self):
        return f"{self.mname}.{self.cname}({self.year}, {self.month}, {self.day})"

    def __format__(self, fmt):
        if not isinstance(fmt, str):
            raise TypeError(f"must be str, not {fmt.__class__.__name__}")
        if len(fmt) != 0:
            return self.strftime(fmt)
        return str(self)

    def __sub__(self, other):
        # if isinstance(other, int):
        #     mjd = self.mjd - other
        #     if mjd in self.valid_mjd:
        #         return Date(*RMJD(mjd))
        #     raise OverflowError
        if isinstance(other, TimeDelta):
            days, fraction = divmod(other, self.resolution)
            if fraction:
                print(f"Warning: {other} has nonzero fraction in date arithmetic")
            mjd = self.mjd - int(days)
            if mjd in self.valid_mjd:
                return type(self).frommjd(mjd)
            raise OverflowError
        if isinstance(other, datetime.date):
            other = Date(other.year, other.month, other.day)
        if isinstance(other, Date):
            return TimeDelta(days=self.mjd - other.mjd)
        return NotImplemented

    def __add__(self, other):
        # if isinstance(other, int):
        #     mjd = self.mjd + other
        #     if mjd in self.valid_mjd:
        #         return Date(*RMJD(mjd))
        #     raise OverflowError
        if isinstance(other, TimeDelta):
            days, fraction = divmod(other, self.resolution)
            if fraction:
                print(f"Warning: {other} has nonzero fraction in date arithmetic")
            mjd = self.mjd + int(days)
            if mjd in self.valid_mjd:
                return type(self).frommjd(mjd)
            raise OverflowError
        return NotImplemented

    __radd__ = __add__

    def __eq__(self, other):
        if isinstance(other, Date) and not isinstance(other, DateTime):
            return self.mjd == other.mjd
        return NotImplemented
    
    def __gt__(self, other):
        if isinstance(other, Date) and not isinstance(other, DateTime):
            return self.mjd > other.mjd
        return NotImplemented

    def __hash__(self):
        if self._hashcode == -1:
            self._hashcode = hash(self._getstate())
        return self._hashcode
    
    def _getstate(self):
        return self.mjd,
        # yhi, ylo = divmod(self.year, 256)
        # return bytes([yhi, ylo, self.month, self.day]),

    # TODO check if this weirdness works for neg years, probably not 
    # def __setstate(self, state):
    #    self.year, self.month, self.day = state

    def __reduce__(self):
        return (self.__class__, self._getstate())

    def __setstate(self, mjd):
        self.year, self.month, self.day = RMJD(mjd)
    #     yhi, ylo, self.month, self.day = string
    #     self.year = yhi * 256 + ylo

Date.min = Date(MINYEAR, 1, 1)
Date.max = Date(MAXYEAR, 12, 31)

class DateTime(Date):
    def __new__(cls, yyyy, mm, dd, hh=0, MM=0, ss=0):
        self = Date.__new__(cls, yyyy, mm, dd)
        self.hh = hh
        self.MM = MM
        self.ss = ss
        return self

    def __eq__(self, other):
        if isinstance(other, (Date, DateTime)):
            return self.mjd == other.mjd
        return NotImplemented
    
    def __gt__(self, other):
        if isinstance(other, (Date, DateTime)):
            return self.mjd > other.mjd
        return NotImplemented

class TZInfo:
    def __init__(self, tzinfo):
        self.tzinfo = tzinfo
        
    def utcoffset(self, dt):
        pass
    
    def dst(self, dt):
        pass
    
    def tzname(self, dt):
        pass
    
    def fromutc(self, dt):
        pass


class Time(metaclass=DTMeta):
    t_pat = "[T| ]?"
    h_pat = "[0-1][0-9]|2[0-4]"
    m_pat = "[0-5][0-9]"
    s_pat = "[0-5][0-9]|60"
    f_pat = "[\\.|,]([0-9]+)"
    thms  = f"({t_pat})({h_pat})(?P<semi>:?)({m_pat})(?P=semi)({s_pat})"
    tz    = "(Z)|(?:([+|-])([0-1][0-9]|2[0-4])(?:(?::?)([0-5][0-9]))?)"
    time  = f"(?:{thms})(?:{f_pat})?({tz})?"
    nano  = Decimal("1000000000")

    __slots__ = 'hour', 'minute', 'second', 'nanosecond', 'tzinfo', '_hashcode', 'fold'

    def __new__(cls, hour=0, minute=0, second=0, nanosecond=0, tzinfo=None,
                *, fold=0):
        
        if minute == 0 and second == 0 and nanosecond == 0 \
            and tzinfo == None and fold == 0:
                if isinstance(hour, int) and hour < 0:
                    self = object.__new__(cls)
                    self.__setstate(hour)
                    self._hashcode = -1
                    return self
        hour, minute, second, nanosecond, fold = _check_time_fields(hour, minute, second, nanosecond, fold)
        _check_tzinfo_arg(tzinfo)
        self = object.__new__(cls)
        self.hour = hour
        self.minute = minute
        self.second = second
        self.nanosecond = nanosecond
        self.tzinfo = tzinfo
        self.fold = fold
        self._hashcode = -1
        return self

    @classmethod
    def fromisoformat(self, time_string):
        m = re.fullmatch(self.time, time_string)
        if m:
            print(m.groups())
            t, hour, semi, minute, second, fraction, \
                tz, z, tzsign, tzh, tzm = m.groups()
            hour = int(hour)
            minute = int(minute)
            second = int(second)
            if fraction:
                l = len(fraction)
                fraction = self.nano * Decimal(fraction)
                nanosecond = int(fraction /  pow(10, l))
            else:
                fraction = 0
            if not tz:
                return Time(hour, minute, second, nanosecond)
            if z == 'Z':
                return Time(hour, minute, second, nanosecond, timezone.utc)
            assert(not z)
            tzsign = 1 - 2 * (tzsign == '-')
            tzh = int(tzh) * tzsign if tzh else 0
            tzm = int(tzm) * tzsign if tzm else 0
            dt = timedelta(hours = tzh, minutes = tzm)
            return Time(hour, minute, second, nanosecond, timezone(dt))
        raise (ValueError(f"{self.cname} Error parsing ISO time {time_string}"))

    @cached_property
    def strf(self):
        mjd_zero = MJD(self.year - 1, 12, 31)
        self.yday = int(self.mjd - mjd_zero)
        weekday_zero = weekday_nr(mjd_zero + MJD0)
        first_sun = (- weekday_zero) % 7
        first_mon = (1 - weekday_zero) % 7
        week_num_U = 1 + (self.yday - first_sun) // 7
        week_num_W = 1 + (self.yday - first_mon) // 7
        strf = {}
        strf["A"] = weekday_str(self.mjd + MJD0)
        strf["a"] = strf["A"][:3]
        strf["w"] = f"{weekday_nr(self.mjd + MJD0)}"
        strf["d"] = f"{self.day:02d}"
        strf["e"] = f"{self.day: 2d}"
        strf["B"] = months[self.month]  # TODO: check values
        strf["b"] = strf["B"][:3]
        strf["m"] = f"{self.month:02d}"
        strf["Y"] = f"{self.year:=-04d}"
        strf["y"] = f"{self.year//100:=-02d}"
        strf["H"] = "00"
        strf["I"] = "00"
        strf["p"] = "AM"
        strf["M"] = "00"
        strf["S"] = "00"
        strf["f"] = "0" * 9
        strf["z"] = ""
        strf["j"] = f"{self.yday:03d}"
        strf["U"] = f"{week_num_U:02d}"
        strf["W"] = f"{week_num_W:02d}"
        strf["z"] = ""
        strf["Z"] = ""
        strf["F"] = f"{self.year}-{self.month}-{self.day}"

        iso_year, iso_week, iso_day = self.isocalendar()
        strf["G"] = f"{iso_year:03d}"
        strf["u"] = f"{iso_day:1d}"
        strf["V"] = f"{iso_week:02d}"
        return strf

    def repl(self, escape):  # replace format % escape by value
        code = escape.group(0)[-1]
        if code in self.strf:
            return self.strf[code]
        elif code == '%':
            return '%'
        # there is a theoretical possibility that locale could cause
        # infinite recursion here.
        elif code == 'c':
            return self.strftime(nl_langinfo(D_T_FMT))
        elif code == 'x':
            return self.strftime(nl_langinfo(D_FMT))
        elif code == 'X':
            return self.strftime(nl_langinfo(T_FMT))
        raise (SyntaxError(f"Encountered invalid % escape ({code})"))

    def strftime(self, fmt):
        r = re.sub("%([a-z|A-Z|%])", self.repl, fmt)
        return (r.format())

    def replace(self, **kwargs):  # TODO check args.
        hour = kwargs["hour"] if "hour" in kwargs else self.hour
        minute = kwargs["minute"] if "minute" in kwargs else self.minute
        second = kwargs["second"] if "second" in kwargs else self.second
        nanosecond = kwargs["nanosecond"] if "nanosecond" in kwargs else self.nanosecond
        tzinfo = kwargs["tzinfo"] if "day" in kwargs else self.tzinfo
        fold = kwargs["fold"] if "fold" in kwargs else self.fold
        return type(self)(hour, minute, second, nanosecond, tzinfo, fold)

    def timetuple(self):
        return struct_time(self.year, self.month, self.day, 0, 0, 0,
                           self.weekday(), self.yday, -1)

    def toordinal(self):
        return self.mjd - self.ORD0

    def isoformat(self, timespec='auto'):
        if self.tzinfo == None:
            tz = ""
        else:
            tz = "TZ"
        if timespec == 'auto':
            if self.nanosecond == 0:
                s = f"{self.hour}"
            n = f".{self.nanosecond:09d}"
        else:
            n = ""
        print(self.tzinfo)
        return f"{self.hour:02d}:{self.minute:02d}:{self.second:02d}{n}"

    def ctime(self):
        return f"{self.strf['a']} {self.strf['b']} {self.day:2d} 00:00:00 {self.year}"

    def __str__(self):
        return self.isoformat()

    def __repr__(self):
        return f"Date({self.year}, {self.month}, {self.day})"

    def __format__(self, fmt):
        print("bar", fmt, self.strftime(fmt))
        if not isinstance(fmt, str):
            raise TypeError(f"must be str, not {type(fmt)}")
        if len(fmt) != 0:
            return self.strftime(fmt)
        return str(self)

    def __sub__(self, b):
        if isinstance(b, int):
            return Date(*RMJD(self.mjd - b))
        elif isinstance(b, Date):
            return TimeDelta(days=self.mjd - b.mjd)
        elif isinstance(b, TimeDelta):
            days, fraction = divmod(b, self.resolution)
            if fraction:
                raise (ValueError(
                    f"{b} has nonzero fraction in date arithmetic"))
            return Date(*RMJD(self.mjd - int(days)))
        raise (TypeError(f"Incompatible types for / : <Date>, {type(b)}"))

    def __add__(self, b):
        if isinstance(b, int):
            return Date(*RMJD(self.mjd + b))
        elif isinstance(b, TimeDelta):
            days, fraction = divmod(b, self.resolution)
            if fraction:
                raise (ValueError(
                    f"{b} has nonzero fraction in date arithmetic"))
            return Date(*RMJD(self.mjd + int(days)))
        raise (TypeError(f"Incompatible types for / : <Date>, {type(b)}"))

    def __bool__(self):
        return True

    def __eq__(self, b):
        return self.year == b.year and self.month == b.month and self.day == b.day


Time.min = Time(0, 0, 0, 0)
Time.max= Time(23, 59, 60, 999999999)
Time.resolution = TimeDelta(nanoseconds=1)


def test_isocalendar():
    week_mondays = [
        ((2003, 12, 22), (2003, 52, 1)),
        ((2003, 12, 29), (2004, 1, 1)),
        ((2004, 1, 5), (2004, 2, 1)),
        ((2009, 12, 21), (2009, 52, 1)),
        ((2009, 12, 28), (2009, 53, 1)),
        ((2010, 1, 4), (2010, 1, 1)),
    ]

    test_cases = []
    for cal_date, iso_date in week_mondays:
        base_date = Date(*cal_date)
        # Adds one test case for every day of the specified weeks
        for i in range(7):
            new_date = base_date + TimeDelta(days=i)
            new_iso = iso_date[0:2] + (iso_date[2] + i,)
            test_cases.append((new_date, new_iso))

    for d, ref in test_cases:
        iso = d.isocalendar()
        date = Date.fromisocalendar(*ref)
        if iso != ref:
            print(d, iso, ref)
        assert (iso == ref)
        if date != d:
            print(d, date, ref)
        assert (d == date)


if __name__ == '__main__':
    td = TimeDelta(hours=1, weeks=1, minutes=15, seconds=10, nanoseconds=1, microseconds=20, milliseconds=100)
    td2 = TimeDelta(days=1, nanoseconds=1)
    td3 = datetime.timedelta(days=1, seconds=2, microseconds=3)
    td4 = TimeDelta(days=1, seconds=0, microseconds=1)
    print(td)
    print(repr(td))
    dt = datetime.datetime.now()
    print(dt, dt+td)
    print(td.total_seconds())
    print(type(td))
    print(td+td2)
    print(td+td3)
    print(td-td2)
    print(td-td3)
    print(td3+td)
    print(td3-td)
    print(-td3)
    print(-td4)
    print(-2*td4)
    print(td4*2.5)
    print(td, td/2.5)
    print(td/td3)
    td5 = TimeDelta(days=1)
    print(td5 < td3)
    print(td3.resolution, td4.resolution)
    print(TimeDelta.resolution)
    print(td//td5)
    date1 = Date(2025, 1, 19)
    date2 = Date(2025, 1, 17)
    print(repr(date1))
    delta = date1 - date2
    print(date1, date2, delta)
    print(date2 + delta, date1 + delta, delta + date1)
    print(TimeDelta(days=27.333) <= TimeDelta(27, hours=24 * 0.333))
    print(TimeDelta(days=-365*(2500+4700), seconds=3.1415))
    #print(-(delta*np.longdouble(2)), type((-(delta/2))))
    date1.min
    print(date1.min, date1.max, date1.resolution)
    # print(date1.strftime("Hello %A (%a) world! (%w)%d %B(%b) %Y-%m-%d %j %U %W%%."))
    # print(Date(2025, 1, 19).isoweekday())
    # print(Date(2003, 12, 29).isocalendar())
    # print(Date(2025, 1, 17).isocalendar())
    # print(Date(2025, 1, 17).weekday())
    # print(Date(2025, 1, 17).ctime())
    # print(time.time())
    # print(MJD(1970, 1, 1).)
    d=Date(-1, 1, 1)
    print(d, d.timetuple())
    print(Date.fromisoformat("2020-W01-1"))
    print(Date.fromisoformat("2004-W01-1"))
    #test_isocalendar()
    print(time.gmtime())
    d = (1, 1, 1)
    print(Date(*d).toordinal())
    print(datetime.date(*d).toordinal())
    print(datetime.datetime.today().__format__("%y"))
    print("{:%y %m}".format(Date(2000, 1, 2)))
    minute = TimeDelta(minutes=1)
    t = TimeDelta(minutes=-2, seconds=30)
    print(repr(divmod(t, minute)))
    print(repr(TimeDelta(seconds=1) * 0.6112295), repr(TimeDelta(nanoseconds=611229500)))
    # print(time.tzname, time.localtime().tm_isdst, time.localtime().tm_zone)
    # print(time.timezone, time.localtime().tm_gmtoff)
    # import zoneinfo
    # # print(zoneinfo.available_timezones())
    # print(dir(zoneinfo.ZoneInfo("Europe/Amsterdam")))
    # print("{:%c}".format(Date(2020, 1, 1)))
    # print(nl_langinfo(D_T_FMT))
    # from decimal import Decimal
    # print(Decimal(3.1415) % 1)

    # print()
    # print(JD(9999, 12, 31) - JD(-4712, 1, 1))
    # print(divmod(TimeDelta(hours=24, minutes=1.2), TimeDelta(hours=1)))
    # print(datetime.timedelta(days=-2, seconds=0, microseconds=0))
    # print(TimeDelta(hours=19)*2, TimeDelta(27) // 2)
    # delta = TimeDelta(days=27, hours=5, microseconds=2500)
    # print(delta, delta.total_days(), float(delta))
    # print(Date(9999, 12, 31) - Date(-4712, 1, 1))
    # print(delta.resolution)
    # date = Date(2000, 1, 19)
    # print(date.replace(year=2025, day=20))
    # print(delta.cname)
    #t = Time.fromisoformat("T000000.123-02")
    #print(t)
    # print(t.resolution)
    #tstart = time.time()
    #print(Date.fromisoformat("2184W531"))
    #tstop=time.time()
    #print(tstop-tstart)
    #print(Time.fromisoformat("T19:27:60+10"))
    #print(Date.fromisoformat("10-05-18"))
    args = 12, 34, 56
    orig = TimeDelta(*args)
    print(date1, date1.mjd)
    d = datetime.date(826, 3, 16)
    D = Date(826,3,12)
    print(Date(*RJD(0)[:3]), Date.min)
    tiny = Date.resolution
    d = Date.min
    print("mjd", d.mjd, d.jd())
    d += tiny
    print(d-tiny)
    print(JD(-4712,1,1))
    print(Date.fromisoformat("2025-W01"))
    i = IsoCalendarDate(1,2,3)
    print(i)
    print(MJD(1945, 11, 12))
    print(Date.strptime("19710510 120001.1", "%Y%m%d %H%M%S.%f"))
    t = Time(1,2,3)
    print(t, type(t.resolution))
    print(Date.fromGregorian(1582,10,14))
    cal = Calendar()
    print(Date.fromisocalendar(2013, 1, 7))