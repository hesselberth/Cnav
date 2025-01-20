#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:13:23 2025

@author: Marcel Hesselberth
"""

from cnav.constants import months, MJD0, SPD
from cnav.dtmath import JD, RJD, MJD, RMJD, weekday_str, weekday_nr
import re
import numpy as np
from functools import total_ordering, cached_property
import time
from collections import namedtuple
import datetime
from locale import nl_langinfo, D_T_FMT, D_FMT, T_FMT
from decimal import Decimal
from enum import Enum


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
PRECISION = 9  # 1 ns
# TODO check consistency of decreased precision
# TODO validate checks

IsoCalendarDate = namedtuple("IsoCalendarDate", ["year", "week", "weekday"])
struct_time = namedtuple("struct_time", ["tm_year", "tm_mon", "tm_mday",
                                         "tm_hour", "tm_min", "tm_sec",
                                         "tm_wday", "tm_yday", "tm_isdst"])

_USEC = 1000  # ns
_MSEC = 1000 * _USEC
_SEC = 1000 * _MSEC
_MIN = 60 * _SEC
_HR = 60 * _MIN
_DAY = 24 * _HR


@total_ordering  # TODO: check efficiency
class TimeDelta(Decimal):

    min = -DELTA_MAX_DAYS * _DAY
    max = DELTA_MAX_DAYS * _DAY

    def __new__(cls, nanoseconds=0, microseconds=0, milliseconds=0,
                seconds=0, minutes=0, hours=0, days=0):
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
        ns = round(ns)
        return super(TimeDelta, cls).__new__(cls, ns)

    def __init__(self, nanoseconds=0, microseconds=0, milliseconds=0,
                 seconds=0, minutes=0, hours=0, days=0):

        days, ns = Decimal.__divmod__(self, _DAY)
        hours, ns = Decimal.__divmod__(ns, _HR)
        minutes, ns = Decimal.__divmod__(ns, _MIN)
        seconds, ns = Decimal.__divmod__(ns, _SEC)
        milliseconds, ns = Decimal.__divmod__(ns, _MSEC)
        microseconds, ns = Decimal.__divmod__(ns, _USEC)
        ns = round(ns)

        self.days = int(days)
        self.hours = int(hours)
        self.minutes = int(minutes)
        self.seconds = int(seconds)
        self.milliseconds = int(milliseconds)
        self.microseconds = int(microseconds)
        self.nanoseconds = int(ns)

    @property
    def resolution(self):
        return TimeDelta(1)

    def total_seconds(self):
        return Decimal.__truediv__(self, _SEC)

    def total_days(self):
        return Decimal.__truediv__(self, _DAY)

    def __str__(self):
        l = []
        ns = 1000000 * self.milliseconds + 1000 * self.microseconds + self.nanoseconds
        if self.days:
            l.append(f"{self.days} days, ")
        l.append(f"{self.hours}:{self.minutes:02d}:{self.seconds:02d}")
        if ns:
            fs = ".{:0" + f"{PRECISION}" + "d}"
            l.append(fs.format(ns))
        return "".join(l)

    def __repr__(self):
        f = []
        if self.days:
            f.append(f"days={self.days}")
        if self.hours:
            f.append(f"hours={self.hours}")
        if self.minutes:
            f.append(f"minutes={self.minutes}")
        if self.seconds:
            f.append(f"seconds={self.seconds}")
        if self.milliseconds:
            f.append(f"milliseconds={self.milliseconds}")
        if self.microseconds:
            f.append(f"microseconds={self.microseconds}")
        if self.nanoseconds:
            f.append(f"nanoseconds={self.nanoseconds}")
        if not f:
            f.append("0")
        # f.append(f"total={self}")
        return f"TimeDelta({', '.join(f)})"

    def __add__(self, b):
        if isinstance(b, TimeDelta):
            return TimeDelta(Decimal.__add__(self, b))
        elif isinstance(b, Date):
            if Decimal.__mod__(self, _DAY) != 0:
                raise (ValueError(f"{b} has nonzero time in date arithmetic"))
            return b + int(Decimal.__floordiv__(self, _DAY))
        raise (TypeError("Incompatible types for + : <TimeDelta>, {type(b)}"))

    def __sub__(self, b):
        if isinstance(b, TimeDelta):
            return TimeDelta(Decimal.__sub__(self, b))
        raise (TypeError(f"Incompatible types for - : <TimeDelta>, {type(b)}"))

    def __neg__(self):
        return TimeDelta(Decimal.__neg__(self))

    def __mul__(self, b):
        if isinstance(b, TimeDelta):
            raise (
                TypeError(f"Incompatible types for * : <TimeDelta>, {type(b)}"))
        if isinstance(b, (int, float, Decimal)):
            return TimeDelta(Decimal.__mul__(self,  b))
        raise (TypeError(f"Incompatible types for * : <TimeDelta>, {type(b)}"))

    def __floordiv__(self, b):
        if isinstance(b, TimeDelta):
            return int(Decimal.__floordiv__(self, b))
        elif isinstance(b, int):
            return TimeDelta(Decimal.__truediv__(self,  b))
        raise (TypeError(f"Incompatible types for // : <TimeDelta>, {type(b)}"))

    def __truediv__(self, b):
        if isinstance(b, TimeDelta):
            return Decimal.__truediv__(self, b)
        elif isinstance(b, (int, float)):
            return TimeDelta(Decimal.__truediv__(self,  b))
        raise (TypeError(f"Incompatible types for / : <TimeDelta>, {type(b)}"))

    def __mod__(self, b):
        if isinstance(b, TimeDelta):
            return TimeDelta(Decimal.__mod__(self, b))
        raise (TypeError(f"Incompatible types for % : <TimeDelta>, {type(b)}"))

    def __divmod__(self, b):
        if isinstance(b, TimeDelta):
            div, mod = Decimal.__divmod__(self, b)
            return int(div), TimeDelta(mod)
        raise (
            TypeError(f"Incompatible types for divmod : <TimeDelta>, {type(b)}"))

    def __float__(self):
        return float(Decimal.__truediv__(self, _DAY))

    def __int__(self):
        return int(float(self))


class Date:
    MIN = (-4712,  1,  1)
    MAX = (9999, 12, 31)
    resolution = TimeDelta(days=1)
    ORD0 = MJD(0, 12, 31)
    dateformats = "([+-]{0,1})(\\d{4})(\\d{2})(\\d{2})|" \
                  "([+|-]{0,1})(\\d{4})-(\\d{2})-(\\d{2})|" \
                  "(?P<s>[+-]{0,1})(\\d{4})-[Ww](\\d{2})-(\\d)"
    datepattern = re.compile(dateformats)

    def __init__(self, year=2000, month=1, day=1):
        self.year = year
        self.month = month
        self.day = day
        self.mjd = MJD(year, month, day)

    @cached_property
    def min(self):
        return Date(*self.MIN)

    @cached_property
    def max(self):
        return Date(*self.MAX)

    @classmethod
    def today(self):
        return self.fromtimestamp(time.time())

    @classmethod
    def fromtimestamp(self, timestamp):
        seconds = np.longdouble(timestamp)
        days = seconds / SPD
        mjd = int(days) + MJD_UNIX_EPOCH
        return Date(*RMJD(mjd))

    @classmethod
    def fromjd(self, jd):
        return Date(*RJD(jd)[:3])

    @classmethod
    def frommjd(self, jd):
        return Date(*RMJD(jd))

    @classmethod
    def fromisoformat(self, date_string):
        m = self.datepattern.match(date_string)
        if m:
            groups = m.groups()
            if groups[1] != None:
                year = int(groups[1]) * (1 - 2 * (groups[0] == '-'))
                month = int(groups[2])
                day = int(groups[3])
                return Date(year, month, day)
            elif groups[5] != None:
                year = int(groups[5]) * (1 - 2 * (groups[4] == '-'))
                month = int(groups[6])
                day = int(groups[7])
                return Date(year, month, day)
            elif groups[9] != None:
                print("g", groups[9])
                year = int(groups[9]) * (1 - 2 * (groups[8] == '-'))
                week = int(groups[10])
                day = int(groups[11])
                return Date.fromisocalendar(year, week, day)
        raise (ValueError(f"Error parsing ISO date {date_string}"))

    @classmethod
    def fromisocalendar(self, year, week, day):
        thursday_year = year
        weekday_1 = Date(year, 1, 1).isoweekday()
        first_thursday = (4 - weekday_1) % 7 + 1
        first_thursday_date = Date(year, 1, first_thursday)
        iso_oldyear = first_thursday_date - 4
        # print(weekday_1, first_thursday, first_thursday_date, iso_oldyear)
        return iso_oldyear + (7 * (week - 1) + day)

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
        strf["f"] = "0" * PRECISION
        strf["z"] = ""
        strf["j"] = f"{self.yday:03d}"
        strf["U"] = f"{week_num_U:02d}"
        strf["W"] = f"{week_num_W:02d}"

        iso_year, iso_week, iso_day = self.isocalendar()
        strf["G"] = f"{iso_year:03d}"
        strf["u"] = f"{iso_day:1d}"
        strf["V"] = f"{iso_week:02d}"
        return strf

    def repl(self, escape):  # replace format % escape by value
        code = escape.group(0)[-1]
        if code in self.strf:
            return (self.strf[code])
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

    def replace(self, **kwargs):
        year = kwargs["year"] if "year" in kwargs else self.year
        month = kwargs["month"] if "month" in kwargs else self.month
        day = kwargs["day"] if "day" in kwargs else self.day
        return Date(year, month, day)

    def timetuple(self):
        return struct_time(self.year, self.month, self.day, 0, 0, 0,
                           self.weekday(), self.yday, -1)

    def toordinal(self):
        return self.mjd - self.ORD0

    def weekday(self):
        return (weekday_nr(self.mjd + MJD0) - 1) % 7

    def isoweekday(self):
        return (weekday_nr(self.mjd + MJD0) - 1) % 7 + 1

    def isocalendar(self):
        iso_weekday = (weekday_nr(self.mjd + MJD0) - 1) % 7 + 1
        thursday_date = self + (4 - iso_weekday)
        iso_year = thursday_date.year
        weekday_1 = Date(iso_year, 1, 1).isoweekday()
        first_thursday = (4 - weekday_1) % 7 + 1
        first_thursday_date = Date(iso_year, 1, first_thursday)
        iso_weeknum = (thursday_date - first_thursday_date).days // 7 + 1
        return IsoCalendarDate(iso_year, iso_weeknum, iso_weekday)

    def isoformat(self):
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"

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

    def __eq__(self, b):
        return self.year == b.year and self.month == b.month and self.day == b.day


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
    date1 = Date(2025, 1, 19)
    date2 = Date(2025, 1, 17)
    print(repr(date1))
    delta = date1 - date2
    print(delta)
    print(date2 + delta, date1 + delta, delta + date1)
    print(TimeDelta(days=27.333) <= TimeDelta(27, hours=24 * 0.333))
    print(TimeDelta(days=-365*(2500+4700), seconds=3.1415))
    # print(-(delta*np.longdouble(2)), type((-(delta/2))))
    print(date1.min, date1.max, date1.resolution)
    print(date1.strftime("Hello %A (%a) world! (%w)%d %B(%b) %Y-%m-%d %j %U %W%%."))
    print(Date(2025, 1, 19).isoweekday())
    print(Date(2003, 12, 29).isocalendar())
    print(Date(2025, 1, 17).isocalendar())
    print(Date(2025, 1, 17).weekday())
    print(Date(2025, 1, 17).ctime())
    print(time.time())
    print(MJD(1970, 1, 1))
    print(Date.today())
    print(Date.fromisoformat("2020-W01-1"))
    print(Date.fromisoformat("2004-w01-1"))
    test_isocalendar()
    print(time.gmtime())
    d = (1, 1, 1)
    print(Date(*d).toordinal())
    print(datetime.date(*d).toordinal())
    print(datetime.datetime.today().__format__("%y"))
    print("{:%y %m}".format(Date(-2000, 1, 2)))
    print(time.tzname, time.localtime().tm_isdst, time.localtime().tm_zone)
    print(time.timezone, time.localtime().tm_gmtoff)
    import zoneinfo
    # print(zoneinfo.available_timezones())
    print(dir(zoneinfo.ZoneInfo("Europe/Amsterdam")))
    print("{:%c}".format(Date(2020, 1, 1)))
    print(nl_langinfo(D_T_FMT))
    from decimal import Decimal
    print(Decimal(3.1415) % 1)

    print()
    print(JD(9999, 12, 31) - JD(-4712, 1, 1))
    print(divmod(TimeDelta(hours=24, minutes=1.2), TimeDelta(hours=1)))
    print(datetime.timedelta(days=-2, seconds=0, microseconds=0))
    print(TimeDelta(hours=19)*2, TimeDelta(27) // 2)
    delta = TimeDelta(days=27, hours=5, microseconds=2500)
    print(delta, delta.total_days(), float(delta))
    print(Date(9999, 12, 31) - Date(-4712, 1, 1))
    print(delta.resolution)
    date = Date(2000, 1, 19)
    print(date.replace(year=2025, day=20))
