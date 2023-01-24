""" RHEAS module for retrieving rainfall data from the Climate Hazard Group
    InfraRed Precipitation with Station (CHIRPS) data archive.

.. module:: chirps
   :synopsis: Retrieve CHIRPS-GEFS rainfall data

.. moduleauthor:: Kostas Andreadis <kandread@jpl.nasa.gov>

"""

from datetime import timedelta
import sys
from . import datasets
from .decorators import geotiff, http
from datetime import datetime


table = "precip.chirpsgefs"


@geotiff
@http
def fetch(dbname, dts, dt, bbox):
    """Downloads CHIRPS-GEFS rainfall forecast data from the data server."""
    """assumes you only want to download a 16-day forecast that was generated at the start date"""
    start = dts[0]
    startday = start.day
    startmonth = start.month
    startyr = start.year
    day = dt.day
    month = dt.month
    yr = dt.year
    url = "https://data.chc.ucsb.edu/products/CHIRPS-GEFS_precip_v12/daily_16day/{0}/{1}/{2}/data.{3}.{4}{5}.tif".format(str(startyr),str(startmonth).zfill(2),str(startday).zfill(2),str(yr),str(month).zfill(2),str(day).zfill(2))
    print(url)
    #url = "ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/p05/{0:04d}/chirps-v2.0.{0:04d}.{1:02d}.{2:02d}.tif.gz"
    return url, bbox, dt



def download(dbname, dts, bbox=None):
    res = 0.05
    for dt in [dts[0] + timedelta(tt) for tt in range((dts[-1] - dts[0]).days + 1)]:
        data, lat, lon, t = fetch(dbname, dts, dt, bbox)
        datasets.ingest(dbname, table, data, lat, lon, res, t)


def dates(dbname):
    dts = datasets.dates(dbname, table)
    return dts
