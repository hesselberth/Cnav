#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 20:28:19 2023

@author: hessel
"""

import os
from re import compile
from time import time
from datetime import datetime
from cnav.dtmath import JD, MJD, RJD, RMJD, bisect
from configparser import ConfigParser


# TODO return zero when dut or leaps not available

path, ext = os.path.splitext(__file__)
config_filename = f"{path}.ini"
config = ConfigParser()
config.read(config_filename)

# Finals
finals_URL = config["Finals"]["URL"]
finals_age = config["Finals"]["maxage"]
finals_dir = config["Finals"]["dir"]

# Leap seconds
leap_URL    = config["Leap seconds"]["URL"]
leap_age    = config["Leap seconds"]["maxage"]
leap_dir    = config["Leap seconds"]["dir"]


class URL_data:
    def __init__(self, url, path, maxage=0, binary=False):
        self.url = url
        self.filename = url.rsplit("/", 1)[-1]
        self.path = path
        self.file_path = os.path.join(self.path, self.filename)
        self.max_file_age = float(maxage) * 86400  # maxage is given in days
        if binary:
            self.readmode  = "rb"
            self.writemode = "wb"
        else:
            self.readmode  = "r"
            self.writemode = "w"
        self.initialize()


    def initialize(self):
        if not os.path.isdir(self.path):
            print(f"Creating path {self.path}.")
            os.makedirs(path)
        if not os.path.isfile(self.file_path):
            self.download()
        if os.path.isfile(self.file_path):
            mtime = os.stat(self.file_path).st_mtime
            file_age = time() - mtime
            if file_age < 0:
                print(f"File age not valid. Please check system clock.")
            else:
                if file_age > self.max_file_age and self.max_file_age > 0:
                    print(f"{self.filename} is out of date.")
                    self.download()
        try:
            file = open(self.file_path, self.readmode)
            data = file.read()
            file.close()
        except:
            raise(FileNotFoundError(f"{self.filename}"))
        else:
            self.decode(data)


    def download(self):
        print(f"Downloading {self.url} to {self.path} ...")
        if os.path.isfile(self.file_path):
            dt = datetime.now().strftime("%Y%m%d%H%M%S")
            directory, filename = os.path.split(self.file_path)
            name, ext = os.path.splitext(filename)
            backup_filename = f"{filename}-{dt}.{ext}"
            backup_path = os.path.join(directory, backup_filename)
            os.rename(self.file_path, backup_path)
            print("Backed up {self.file_path} to {backup_path} .")
        import urllib.request
        urllib.request.urlretrieve(self.url, self.file_path)


    def decode(self, txt):
        pass


    def __str__(self):
        return self.file_path


class Finals(URL_data):
    def __init__(self, url, path=".", backup_path="backup", maxage=0):
        URL_data.__init__(self, url, path, backup_path, maxage)
        print(self)

    def decode(self, txt):
        self.info = {}
        self.mjd = []
        self.prev_mjd = 0
        self.pred = False
        self.filedate = None
        lines = txt.split('\n')
        for line in lines:
            try:
                yy    = int(line[0:2])
                mm    = int(line[2:4])
                dd    = int(line[4:6])
                mjd   = float(line[7:15])
                pred  = line[57].upper()
                dut1  = float(line[58:68])
                error = float(line[68:78])
            except:
                pass
            else:
                self.add_dut1(yy, mm, dd, mjd, pred, dut1, error)

    def add_dut1(self, yy, mm, dd, mjd, pred, dut1, error):
        if pred == 'P':
            pred = True
        elif pred == 'I':
            pred = False
        else:
            raise(SyntaxError("Entry is neither IERS data nor prediction."))
        if mjd >= 51544:        # Fix date from a 2-digit year
            yyyy = yy + 2000
        else:
            yyyy = yy + 1900
        assert(mjd == MJD(yyyy, mm, dd))
        if mjd > self.prev_mjd:
            self.prev_mjd = mjd
        else:
            raise(ValueError("Dates in IAU 2000 file do not monotonically increase."))
        self.mjd.append(mjd)
        self.info[mjd] = dut1, error
        if pred and not self.pred:  # first prediction, save file date
            self.pred = True
            self.filedate = mjd - 1

    def __str__(self):
        if not self.filedate:
            return ""
        today = datetime.now()
        mjd_today = MJD(today.year, today.month, today.day)
        try:
            y0, m0, d0, t = RMJD(self.filedate)
            dut1, error = self.info[mjd_today]
            ymax, mmax, dmax, t = RMJD(self.mjd[-1])
        except:
            return "Failed to retrieve DUT1 data"
        errstr = "\nDUT1 error today (%d-%0.2d-%0.2d) is %.3f seconds." % \
                (today.year, today.month, today.day, error)
        return "DUT1 data from %d-%0.2d-%0.2d. Data available until %d-%0.2d-%0.2d.%s" % \
                (y0, m0, d0, ymax, mmax, dmax, errstr)

    def __call__(self, mjd):  # DUT1 in seconds, mjd from UTC 
        t = mjd % 1
        mjd = mjd // 1
        try:
            dut1_prev, error = self.info[mjd]
            dut1_next, error = self.info[mjd+1]
        except KeyError:
            print("Warning: DUT1 value not in table. Using 0.")
            return 0
        else:
            if dut1_next - dut1_prev > 0.5:  # next day has positive leap second
                dut1_next -= 1.0
            assert(dut1_next > -1)
            if t > 1.0:
                t = 1.0
            return dut1_prev + (dut1_next - dut1_prev) * t


class Leapseconds(URL_data):
    fmt =  "([0-9]{4})\\s+([A-Z]{3})[\\s]+([0-9]{1,2})[\\s]+=JD\\s([0-9.]+)"
    fmt += "[\\s]+TAI-UTC=[\\s]*([0-9.]+)[\\s]*(S)[\\s]*\\+[\\s]*"
    fmt += "\\([\\s]*MJD[\\s]*-[\\s]*([0-9.]+)[\\s]*\\)[\\s]*X[\\s]*([0-9.]+)"
    fmt += "[\\s]*(S)[\\s]*\n"
    months = { "jan":1, "feb":2, "mar":3, "apr":4, "may":5, "jun":6,   \
              "jul":7, "aug":8, "sep":9, "oct":10, "nov":11, "dec":12 }

    def __init__(self, url, path=".", backup_path="backup", maxage=0):
        URL_data.__init__(self, url, path, backup_path, maxage)
        print(self)

    def decode(self, txt):
        self.info = {}
        self.jd = []
        self.rinfo = {}
        self.rjd = []
        self.prevjd = 0
        pattern = compile(self.fmt)
        leaps = pattern.findall(txt)
        for leap in leaps:
            y, m, d, jd, offset, s1, mjd, factor, s2 = leap
            self.add_leap(y, m, d, jd, offset, s1, mjd, factor, s2)

    def add_leap(self, y, m, d, jd, offset, s1, mjd, factor, s2):
        if s1.upper() != 'S' or s2.upper() != 'S':
            raise(SyntaxError("unknown unit in leap second file "
                              "(seconds expected)."))
        yyyy = int(y)
        mm   = self.months[m.lower()]
        dd   = int(d)
        jd   = float(jd)
        assert(jd == JD(yyyy, mm, dd))
        offset  = float(offset)
        mjd     = float(mjd)
        factor  = float(factor)
        if jd > self.prevjd:
            self.prevjd = jd
        else:
            raise(ValueError("Dates in leap second file do not "
                             "monotonically increase."))
        self.jd.append(jd)
        self.info[jd] = offset, mjd, factor
        yr, mr, dr, tr = RJD(jd)

        assert(yyyy == yr and mm == mr and dd == dr)
        # checks invariance of reversejd(jd())

    def __str__(self):
        lastleap = self.jd[-1]
        yyyy, mm, dd, t = RJD(lastleap)      # to calculate value
        value = self(yyyy, mm, dd)
        yyyy, mm, dd, t = RJD(lastleap - 1)  # to calculate date
        return "Last leap second was at %d-%0.2d-%0.2d. " \
               "TAI-UTC = %.3f seconds." % (yyyy, mm, dd, value)

    def new_leap_second(self, YYYY, MM, DD):  # day starts with new leap second
        jd = JD(YYYY, MM, DD)
        return (jd in self.jd)

    def __call__(self, YYYY, MM, DD):  # TAI - UTC in seconds, date UTC
        #F = (hh * 3600 + mm * 8 + ss) / 86400
        jd = JD(YYYY, MM, DD)
        mjd = jd - 2400000.5
        last_leap = bisect(self.jd, jd)  # includes today
        if last_leap:
            offset, mjd0, factor = self.info[last_leap]
            return offset + (mjd - mjd0) * factor
        return 0

def finals():
    return Finals(finals_URL, finals_dir, finals_age)

def leapseconds():
    return Leapseconds(leap_URL, leap_dir, leap_age)
