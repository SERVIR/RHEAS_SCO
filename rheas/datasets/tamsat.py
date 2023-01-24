""" RHEAS module for retrieving rainfall data from the TAMSAT data archive.

.. module:: tamsat
   :synopsis: Retrieve TAMSAT rainfall data

.. moduleauthor:: Sara Miller <sem0029@uah.edu>
adapted from CHIRPS data download script by Kostas Andreadis <kandread@umass.edu>

"""

from datetime import timedelta

from . import datasets
from .decorators import netcdftam, http

table = "precip.tamsat"


@netcdftam
@http
def fetch(dbname, dt, bbox):
    """Downloads TAMSAT rainfall data from the data server."""
    url = "https://gws-access.jasmin.ac.uk/public/tamsat/rfe/data/v3.1/daily/{0:04d}/{1:02d}/rfe{0:04d}_{1:02d}_{2:02d}.v3.1.nc"
    return url, bbox, dt


def download(dbname, dts, bbox=None):
    res = 0.05
    for dt in [dts[0] + timedelta(tt) for tt in range((dts[-1] - dts[0]).days + 1)]:
        data, lat, lon, t = fetch(dbname, dt, bbox)
        t= dt
        datasets.ingest(dbname, table, data, lat, lon, res, t)


def dates(dbname):
    dts = datasets.dates(dbname, table)
    return dts
