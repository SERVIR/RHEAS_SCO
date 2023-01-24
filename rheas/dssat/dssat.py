""" Class definition for the DSSAT model interface

.. module:: dssat
   :synopsis: Definition of the DSSAT model class

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>
edited by Narendra Das to include fertilizers and irrigation from the database
edited by Sara Miller to remove the planting dates read from database and take less time to pull fertilizers from db

"""

import decimal
import distutils.core
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import date, datetime, timedelta

import numpy as np
import pandas as pd

from .. import dbio

class DSSAT(object):
    def __init__(self, dbname, name, resolution, startyear, startmonth, startday, endyear, endmonth, endday, nens, vicopts, shapefile=None,fertilizer=None,irrigation=None,plantingdate=None, assimilate=True):
        log = logging.getLogger(__name__)
        #self.path = tempfile.mkdtemp(dir=".")
        self.path = 'tmpD'+datetime.now().strftime('%Y%m%d%H%M%S')
        if not os.path.isdir(self.path):
            os.mkdir(self.path)
        self.startyear = startyear
        self.startmonth = startmonth
        self.startday = startday
        self.endyear = endyear
        self.endmonth = endmonth
        self.endday = endday
        self.crop = None
        self.cultivars = {}
        self.lat = []
        self.lon = []
        self.elev = []
        self.depths = []
        self.dbname = dbname
        self.name = name
        self.res = resolution
        self.nens = nens
        self.shapefile = shapefile
        self.assimilate = assimilate
        self.plantingdate=plantingdate
        self.fertilizer=fertilizer
        self.irrigation=irrigation
        self.data_path = os.path.dirname(os.path.abspath(__file__)) + "/../../data"
        if not os.path.isdir(self.data_path):
            try:
                self.data_path = os.environ['RHEASDATA']
            except KeyError:
                log.error("Data directory is not accessible. Please set RHEASDATA environment variable. Exiting...")
                sys.exit()
        self.modelpaths = {}
        self.modelstart = {}
        self.grid_decimal = - (decimal.Decimal(str(self.res)).as_tuple().exponent - 1)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        if 'lai' in vicopts or ('save' in vicopts and vicopts['save'].find("lai") >= 0):
            self.lai = "vic"
        else:
            self.lai = None
        if 'save to' in vicopts:
            self.datafrom = vicopts['save to']
        else:
            self.datafrom = "db"
        cur.execute(
            "select * from information_schema.tables where table_name='basin' and table_schema=%s", (name,))
        if not bool(cur.rowcount):
            log.error("No simulation named {0} exists in database. You might have to run VIC.".format(name))
            sys.exit()
        cur.execute(
            'select basefile from vic.input where resolution=%f;' % self.res)
        self.basefile = "{0}/{1}".format(self.data_path, cur.fetchone()[0])
        cur.close()
        db.close()

    def readVICSoil(self):
        """Extract information from VIC database table on latitude, longitude,
        elevation  and soil depths."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sql = "select st_y(geom), st_x(geom), elev, depths from {0}.basin".format(
            self.name)
        cur.execute(sql)
        pixels = cur.fetchall()
        self.lat, self.lon, self.elev, self.depths = zip(*pixels)
        self.lat = np.array(self.lat)
        self.lon = np.array(self.lon)
        self.elev = np.array(self.elev)
        self.depths = list(self.depths)
        cur.close()
        db.close()

    def writeWeatherFiles(self, modelpath, name, year, month, day, weather, elev, lat, lon, ts=None, te=None):
        """Writes ensemble weather files for specific pixel."""
        if isinstance(weather, list):
            data = (weather * (int(self.nens / len(weather)) + 1))[:self.nens]
        else:
            data = [weather] * self.nens
        for ens in range(self.nens):
            filename = "{0}/WEATH{1:03d}.WTH".format(modelpath, ens + 1)
            fout = open(filename, 'w')
            fout.write("*WEATHER DATA : {0}\r\n".format(name[:5].upper()))
            fout.write("\r\n")
            fout.write("@ INSI LAT LONG ELEV TAV AMP REFHT WNDHT\r\n")
            tavg = np.mean(data[ens][:, 1:3])
            fout.write("{0:6s} {1:6.2f} {2:6.2f} {3:5.0f} {4:5.1f} {5:5.1f} {6:5.1f} {7:5.1f} \r\n".format(
                name[:5].upper(), lat, lon, elev, tavg, -99.0, -99.0, -99.0))
            fout.write("@DATE SRAD TMAX TMIN RAIN DEWP WIND PAR\r\n")
            if ts is None or te is None:
                ts = 0
                te = len(data[ens])
            for p in range(ts, te):
                datestr = str(int(year[p]))[-2:] + date(int(year[p]),
                                                        int(month[p]), int(day[p])).strftime("%j")
                fout.write("{0}  {1:4.1f}  {2:4.1f}  {3:4.1f}  {4:4.1f}\r\n".format(
                    datestr, data[ens][p, 0] * 0.086400, data[ens][p, 1], data[ens][p, 2], data[ens][p, 3]))
            fout.close()

    def readVICOutputFromFile(self, lat, lon, depths, filespath):
        """Read DSSAT inputs from VIC output files for a specific pixel."""
        startdate = date(self.startyear, self.startmonth, self.startday)
        enddate = date(self.endyear, self.endmonth, self.endday)
        filename = "{0}/output/eb_{1:.{3}f}_{2:.{3}f}".format(
            filespath, lat, lon, self.grid_decimal)
        viceb = np.loadtxt(filename)
        filename = "{0}/output/sub_{1:.{3}f}_{2:.{3}f}".format(
            filespath, lat, lon, self.grid_decimal)
        vicsm = np.loadtxt(filename)
        filename = "{0}/output/sur_{1:.{3}f}_{2:.{3}f}".format(
            filespath, lat, lon, self.grid_decimal)
        vicsr = np.loadtxt(filename)
        filename = "{0}/forcings/data_{1:.{3}f}_{2:.{3}f}".format(
            filespath, lat, lon, self.grid_decimal)
        met = np.loadtxt(filename)
        sm = vicsm[:, 3:len(depths) + 3]
        weather = np.vstack(
            (viceb[:, 3] + viceb[:, 4], met[:, 1], met[:, 2], met[:, 0])).T
        year = vicsm[:, 0].astype(int)
        month = vicsm[:, 1].astype(int)
        day = vicsm[:, 2].astype(int)
        tidx = [i for i in range(len(year)) if date(year[i], month[i], day[
            i]) >= startdate and date(year[i], month[i], day[i]) <= enddate]
        lai = dict(zip([date(year[i], month[i], day[i])
                        for i in range(len(year)) if i in tidx], vicsr[:, 12]))
        return year[tidx], month[tidx], day[tidx], weather[tidx, :], sm[tidx, :], lai

    def readVICOutputFromDB(self, gid, depths, planting, harvest):
        """Read DSSAT inputs from database."""
        if planting is None:
            startdate = date(self.startyear, self.startmonth, self.startday)
        else:
            startdate = planting       
        if harvest is None:
            enddate = date(self.endyear, self.endmonth, self.endday)            
        else:
            enddate = harvest
        ndays = (enddate - startdate).days + 1  
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        date_sql = "fdate>=date '{0}' and fdate<=date '{1}'".format(startdate.strftime("%Y-%m-%d"), enddate.strftime("%Y-%m-%d"))       
        data = {}
        ii=0
        varnames = ["net_short", "net_long",
                    "soil_moist", "rainf", "tmax", "tmin"]
        if self.lai is not None:
            varnames.append("lai")
        else:
            lai = None            
        ii=ii+1
        for varname in varnames:
            sqlvars = ["fdate"]
            sql = "select column_name from information_schema.columns where table_schema='{0}' and table_name='{1}' and column_name='ensemble'".format(
                self.name, varname)
            cur.execute(sql)
            if bool(cur.rowcount):
                sqlvars += ["ensemble"]
            sql = "select column_name from information_schema.columns where table_schema='{0}' and table_name='{1}' and column_name='layer'".format(
                self.name, varname)
            cur.execute(sql)
            if bool(cur.rowcount):
                sqlvars += ["layer"]
            sql = "select max(fdate) from {0}.{1} where {2} and ensemble=0".format(self.name, varname, date_sql)
            cur.execute(sql)
            nowcast_date = cur.fetchall() 
            nowcast_date = nowcast_date[0][0]
            if (nowcast_date >= enddate):
                nowcast_date = enddate
            date_upto_nowcast_sql = "fdate>=date '{0}' and fdate<=date '{1}'".format(startdate.strftime("%Y-%m-%d"), nowcast_date.strftime("%Y-%m-%d"))
            if (nowcast_date + timedelta(1) <= enddate):                
                date_after_nowcast_sql = "fdate>=date '{0}' and fdate<=date '{1}'".format((nowcast_date + timedelta(1)).strftime("%Y-%m-%d"), enddate.strftime("%Y-%m-%d"))      
                sql = "select {0}, avg((st_summarystats(rast)).mean) from {1}.{2}, {1}.agareas where st_intersects(rast,geom) and gid={3} and {4} group by gid,{0} order by fdate".format(",".join(sqlvars), self.name, varname, gid, date_after_nowcast_sql)
                cur.execute(sql)
                results_after_nowcast = cur.fetchall()
                                 
            sql = "select {0}, avg((st_summarystats(rast)).mean) from {1}.{2}, {1}.agareas where st_intersects(rast,geom) and {3} group by {0} order by fdate".format(",".join(sqlvars), self.name, varname, date_upto_nowcast_sql)
            cur.execute(sql)
            results_upto_nowcast = cur.fetchall()
            #print('***********************',results_upto_nowcast)          
            if "ensemble" in sqlvars:
                if (nowcast_date == enddate):                
                    vicnens = np.max([r[1] for r in results_upto_nowcast])                   
                else:               
                    vicnens = np.max([r[1] for r in results_after_nowcast])
                if vicnens == 0:                  
                    data[varname] = np.array([r[-1] for r in results_upto_nowcast])
                    if "layer" in sqlvars:
                        layers = np.array([r[2] for r in results_upto_nowcast])
                        nlayers = np.max(layers)
                    else:
                        year = np.array([r[0].year for r in results_upto_nowcast])        
                        month = np.array([r[0].month for r in results_upto_nowcast])           
                        day = np.array([r[0].day for r in results_upto_nowcast])
                else:        
                    #data[varname] = [np.array([r[-1] for r in results if r[1] == ens + 1 or r[1] == 0]) for ens in range(vicnens)]
                    data_array_ancst = [np.array([r[-1] for r in results_after_nowcast if r[1] == ens + 1]) for ens in range(vicnens)]
                    data_array_uncst = np.array([r[-1] for r in results_upto_nowcast])
                    for ens in range(vicnens):
                        temp_arr = data_array_ancst[ens]
                        temp_arr = np.concatenate((data_array_uncst, temp_arr), axis=0)
                        data_array_ancst[ens] = temp_arr
                    data[varname] = data_array_ancst            
                    if "layer" in sqlvars:                      
                        layers1 = np.array([r[2] for r in results_after_nowcast if r[1] == 1 or r[1] == 0])
                        layers2 = np.array([r[2] for r in results_upto_nowcast if r[1] == 1 or r[1] == 0])
                        layers = np.concatenate((layers2, layers1), axis=0)
                        nlayers1 = np.max(layers1)
                        nlayers2 = np.max(layers2)
                        nlayers = max(nlayers1, nlayers2)
                    else:                
                        year1 = np.array([r[0].year for r in results_after_nowcast if r[1] == 1 or r[1] == 0])
                        month1 = np.array([r[0].month for r in results_after_nowcast if r[1] == 1 or r[1] == 0])
                        day1 = np.array([r[0].day for r in results_after_nowcast if r[1] == 1 or r[1] == 0])
                        year2 = np.array([r[0].year for r in results_upto_nowcast if r[1] == 1 or r[1] == 0])
                        month2 = np.array([r[0].month for r in results_upto_nowcast if r[1] == 1 or r[1] == 0])
                        day2 = np.array([r[0].day for r in results_upto_nowcast if r[1] == 1 or r[1] == 0])
                        year = np.concatenate((year2, year1), axis=0)
                        month = np.concatenate((month2, month1), axis=0)
                        day = np.concatenate((day2, day1), axis=0)
            else:       
                data[varname] = np.array([r[-1] for r in results_upto_nowcast])
                if "layer" in sqlvars:
                    layers = np.array([r[1] for r in results_upto_nowcast])
                    nlayers = np.max(layers)
                else:
                    year = np.array([r[0].year for r in results_upto_nowcast])
                    month = np.array([r[0].month for r in results_upto_nowcast])
                    day = np.array([r[0].day for r in results_upto_nowcast])
                # assert len(year) == ndays and len(month) == ndays and len(day) == ndays
        cur.close()
        db.close()
        if "ensemble" in sqlvars:           
            if vicnens == 0:
                weather = np.vstack(
                    (data["net_short"] - data["net_long"], data["tmax"], data["tmin"], data["rainf"])).T            
                if self.lai is not None:
                    lai = dict(zip([date(year[i], month[i], day[i])
                                    for i in range(len(year))], np.array(data["lai"]).T))
                sm = np.zeros((len(year), nlayers))
                for l in range(nlayers):
                    sm[:, l] = [m for mi, m in enumerate(
                        data["soil_moist"]) if layers[mi] == l + 1]
            else:
                weather = [np.vstack((data["net_short"][e] - data["net_long"][e], data["tmax"][
                              e], data["tmin"][e], data["rainf"][e])).T for e in range(len(data["net_short"]))]
                if self.lai is not None:
                    lai = dict(zip([date(year[i], month[i], day[i]) for i in range(
                        len(year))], np.mean(np.array(data["lai"]).T, axis=1)))   
                sm = [np.zeros((len(year), nlayers))] * len(data["soil_moist"])
                for e in range(len(sm)):
                    for l in range(nlayers):
                        sm[e][:, l] = [m for mi, m in enumerate(
                            data["soil_moist"][e]) if layers[mi] == l + 1]
        else:
            weather = np.vstack(
                (data["net_short"] - data["net_long"], data["tmax"], data["tmin"], data["rainf"])).T
            if self.lai is not None:
                lai = dict(zip([date(year[i], month[i], day[i])
                                for i in range(len(year))], np.array(data["lai"]).T))
            sm = np.zeros((len(year), nlayers))
            for l in range(nlayers):
                sm[:, l] = [m for mi, m in enumerate(
                    data["soil_moist"]) if layers[mi] == l + 1]
        return year, month, day, weather, sm, lai

    def readVICOutput(self, gid, depths, planting=None, harvest=None):
        """Reads DSSAT time-varying inputs by reading either from files or a database."""
        log = logging.getLogger(__name__)
        if isinstance(self.datafrom, list):
            inputs = []
            while len(inputs) < self.nens:
                inputs += self.datafrom
            inputs = inputs[:self.nens]
            lat, lon = self.gid[gid]
        if self.datafrom == 'db':
            year, month, day, weather, sm, lai = self.readVICOutputFromDB(
                gid, depths, planting, harvest)  
        else:
            log.error("VIC output was not saved in the database. Cannot proceed with the DSSAT simulation.")
            sys.exit()
        return year, month, day, weather, sm, lai

    def writeLAI(self, modelpath, gid, year, month, day, ts=None, te=None, viclai=None, tablename="lai.modis"):
        """Writes LAI file for DSSAT."""
        fout = open("{0}/LAI.txt".format(modelpath), 'w')
        if ts is None or te is None:
            ts = 0
            te = len(year)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute("select * from information_schema.tables where table_name=%s and table_schema='lai'",
                    (tablename.split(".")[1],))
        if bool(cur.rowcount) and not self.lai == "vic":
            dt1 = date(year[ts], month[ts], day[ts])
            dt2 = date(year[te-1], month[te-1], day[te-1])
            sql = "select fdate,avg((st_summarystats(st_clip(rast,geom))).mean) from {0},{1}.agareas where st_intersects(rast,geom) and fdate>=date '{2}' and fdate<=date '{3}' and gid={4} group by fdate".format(tablename, self.name, dt1.strftime("%Y-%m-%d"), dt2.strftime("%Y-%m-%d"), gid)
            cur.execute(sql)
            if bool(cur.rowcount):
                results = cur.fetchall()
                lai = {}
                for r in results:
                    if r[1] is None:
                        lai[r[0]] = -9999.0
                    else:
                        lai[r[0]] = r[1] / 10.0
            else:
                lai = {}
        else:
            lai = viclai
        for t in range(ts, te):
            dt = date(year[t], month[t], day[t])
            if lai is not None and dt in lai:
                fout.write("{0:.1f}\n".format(lai[dt]))
            else:
                fout.write("-9999.0\n")
        fout.close()
        cur.close()
        db.close()

    def writeSoilMoist(self, modelpath, year, month, day, smi, dz, ts=None, te=None):
        """Writes soil moisture information file."""
        filename = "{0}/SOIL_MOISTURE.ASC".format(modelpath)
        fout = open(filename, 'w')
        if ts is None or te is None:
            ts = 0
            te = len(smi)
        for t in range(ts, te):
            dt = date(year[t], month[t], day[t])
            doy = int(dt.strftime("%j"))
            fout.write("{0:.0f} {1:.0f} {2:.0f} ".format(
                dt.year, dt.month, dt.day))
            for lyr in range(len(dz)):
                fout.write("{0:.3f} ".format(smi[t, lyr]))
            fout.write("{0}\n".format(doy))
        fout.close()

    def sampleSoilProfiles(self, gid):
        """Samples soil profiles from database to be used in DSSAT control file."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sql = "with f as (select st_envelope(geom) as geom from {0}.agareas where gid={1}) select props from dssat.soils as s,f where st_intersects(s.geom,f.geom)".format(self.name, gid)
        cur.execute(sql)
        # if crop area is too small, look for nearest soil profiles
        dist = 0.1
        while not bool(cur.rowcount):
            sql = "with a as (select st_buffer(geom,{2}) as geom from {0}.agareas where gid={1}) select props from dssat.soils as s,a where st_intersects(s.geom,a.geom)".format(
                self.name, gid, dist)
            dist += 0.1
            cur.execute(sql)
        profiles = cur.fetchall()
        ens = np.random.choice(range(len(profiles)), self.nens)
        cur.close()
        db.close()
        return [profiles[e] for e in ens]

    def writeConfigFile(self, modelpath, nlayers, startdate, enddate):
        """Write DSSAT-ENKF config file."""
        configfilename = "ENKF_CONFIG.TXT"
        fout = open("{0}/{1}".format(modelpath, configfilename), 'w')
        fout.write("!Start_DOY_of_Simulation:\n{0}\n".format(
            int(startdate.strftime("%j"))))
        fout.write("!End_DOY_of_Simulation\n{0}\n".format(
            int(enddate.strftime("%j"))-1))
        fout.write("!Year_of_Simulation:\n{0}\n".format(startdate.year))
        fout.write("!Ensemble_members\n{0}\n".format(self.nens))
        fout.write("!Number_of_soil_layers\n{0}\n".format(nlayers))
        ndays = (enddate - startdate).days
        fout.write("!Number_of_RS_data\n{0}".format(ndays))
        fout.close()
        return configfilename

    def calcCroplandFract(self):
        """Calculate fraction of cropland for specific pixel."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sql = "select gid,avg((st_summarystats(st_clip(rast,geom))).mean) from dssat.cropland,{0}.agareas where st_intersects(rast,geom) group by gid order by gid".format(
            self.name)
        cur.execute(sql)
        fract = dict((r[0], r[1]) for r in cur.fetchall())
        cur.close()
        db.close()
        return fract

    def readShapefile(self):
        """Read areas from shapefile where DSSAT will be run."""
        log = logging.getLogger(__name__)
        try:
            cmd = "shp2pgsql -s 4326 -d -I -g geom {0} {1}.agareas | psql -d {2}".format(self.shapefile, self.name, self.dbname)
            subprocess.call(cmd, shell=True)
            db = dbio.connect(self.dbname)
            cur = db.cursor()
            sql = "select gid, st_x(st_centroid(geom)), st_y(st_centroid(geom)) from {0}.agareas".format(self.name)
            cur.execute(sql)
            geoms = cur.fetchall()
            return geoms
        except IOError:
            log.error("Shapefile {0} for DSSAT simulation does not exist. Exiting...".format(
                self.shapefile))
            sys.exit()

    def readPlanting(self, lat, lon, fromShapefile=False):
        """Retrieve planting dates for pixel."""
        if self.crop is None:
            self.crop = "maize"
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sql = "select st_nearestvalue(rast,st_geomfromtext('POINT({0} {1})',4326)) as doy from crops.plantstart where type like '{2}' and st_intersects(rast,st_geomfromtext('POINT({0} {1})',4326)) order by doy".format(
            lon, lat, self.crop)
        cur.execute(sql)
        results = cur.fetchall()
        # if this is a multi-year simulation append planting dates for each successive year
        plantdates = []
        for yr in range(self.startyear, self.endyear + 1):
            for r in results:
                if r[0] is not None:
                    plantdates.append(date(yr, 1, 1) + timedelta(r[0] - 1))
        cur.close()
        db.close()
        startdt = date(self.startyear, self.startmonth, self.startday)
        planting = [p for p in plantdates if p >= startdt and p <= date(self.endyear, self.endmonth, self.endday)]
        # with the exception of no planting dates satisfying that condition (e.g. short-term forecasts)
        # in that case, we will find the closest past date to the simulation start date
        aplant = []
        sdt = int(startdt.strftime('%j'))
        for r in results:
            if r[0] is not None:
                if r[0] > sdt:
                    aplant.append(date(startdt.year-1, 1, 1) + timedelta(r[0] - 1))
                else:
                    aplant.append(date(startdt.year, 1, 1) + timedelta(r[0] - 1))
        planting = [max(aplant)] + planting
        #planting=np.full([self.nens,len(planting)],planting)
        
        return planting
   
    def planting(self, lat, lon,fromShapefile=False):
        """Retrieve planting dates for pixel."""
        log = logging.getLogger(__name__)
        if self.crop is None:
            self.crop = "maize"
        db = dbio.connect(self.dbname)
        cur = db.cursor()    
        self.plantingdate='static'
        if not self.plantingdate or self.plantingdate == 'static':          
            planting=self.readPlanting(lat, lon, False)
        
        elif self.plantingdate == 'dynamic':
            #extracting planting date for the domain, output would be the planting date and the number of each date located in the domain
            if dbio.tableExists(self.dbname, "crops", "plantstart"):
                sql = "select ST_ValueCount(rast) from crops.plantstart where fdate like '{2}' and source like'{3}' and type like '{4}' and st_intersects(rast,st_geomfromtext('POINT({0} {1})',4326))".format(
                lon, lat, self.startyear,'sar',self.crop)  
                cur.execute(sql)
                if bool(cur.rowcount):
                    results = cur.fetchall() 
                    pixelcount=[]
                    PD=[]
                    pr=[]    
                    for item in results:
                        x=item[0].split(',')
                        PD.append(int(x[0][1:]))
                        pixelcount.append(int(x[1][:-1]))              
                    PD=np.array(PD)
                    PD_id=np.argsort(PD)
                    pixelcount=np.array(pixelcount)               
                    march_count=np.where((PD>59)&(PD<91),pixelcount,0)
                    april_count=np.where((PD>91)&(PD<121),pixelcount,0)
                    may_count=np.where((PD>121)&(PD<152),pixelcount,0)
                    june_count=np.where((PD>153)&(PD<182),pixelcount,0)
                    july_count=np.where((PD>182)&(PD<213),pixelcount,0)
                    aug_count=np.where((PD>213)&(PD<243),pixelcount,0)          
                    m_PD=[march_count,april_count,may_count,june_count,july_count,aug_count]          
                    pisum=[np.sum(march_count),np.sum(april_count),np.sum(may_count),np.sum(june_count),np.sum(july_count),np.sum(aug_count)]
                    for i in range(len(pisum)):
                        pr.append(np.round(self.nens*pisum[i]/np.sum(pisum),decimals=0))
                    pr=np.array(pr)  
                    a=[]
                    for i in range(len(pr)):
                        for j in range(len(m_PD[0])):
                            if m_PD[i][j]>0:                           
                                prr=(np.round(pr[i]*m_PD[i][j]/pisum[i],decimals=0))
                                a=a+list(np.full(int(prr),int(PD[j]))) 
                    if self.nens-len(a)>0:    
                        a=a+list(np.full(int(self.nens-len(a)),PD[PD_id[0]]))
                    elif self.nens-len(a)<0:
                        a=a[0:self.nens]       
                    # if this is a multi-year simulation append planting dates for each successive year
                    plantdates=[]
                    for yr in range(self.startyear, self.endyear + 1):
                        for r in a:
                            plantdates.append(date(yr, 1, 1) + timedelta(int(r) - 1))
                    cur.close()
                    db.close()
                    startdt = date(self.startyear, self.startmonth, self.startday)
                    planting=plantdates
                #planting = [p for p in plantdates if p >= startdt and p <= date(self.endyear, self.endmonth, self.endday)]
                else:
                    log.warning("Static planting is used because there is no data in the table for this region")
                    planting=self.readPlanting(lat, lon, False)
            else:
                log.warning("Static planting is used because the table does not exist")
                planting=self.readPlanting(lat, lon, False)    
        elif self.plantingdate.isdigit:
            print(self.plantingdate)
            a=np.full(self.nens,int(self.plantingdate))
            plantdates=[]
            for yr in range(self.startyear, self.endyear + 1):
                for r in a:
                    plantdates.append(date(yr, 1, 1) + timedelta(int(r) - 1))
            startdt = date(self.startyear, self.startmonth, self.startday)
            planting=plantdates
            
        else:
            log.warning("Static planting is used because the planting method was not selected correctly")
            planting=self.readPlanting(lat, lon, False)
        return planting
            


    def interpolateSoilMoist(self, sm, depths, dz):
        """Estimate soil moisture at DSSAT depths."""
        sm_i = []
        if len(sm.shape) < 2:
            sm = np.reshape(sm, (1, len(sm)))
        for t in range(sm.shape[0]):
            u = sm[t, :] / np.array(depths * 1000.0)
            z = [100.0 * depths[0] / 2.0]
            for lyr in range(1, len(u)):
                # midpoint of each layer in cm
                z.append(100.0 * (depths[lyr - 1] + depths[lyr] / 2.0))
            dz1 = [0.0] + list(dz)
            znew = np.array([dz1[i] + (dz1[i + 1] - dz1[i]) /
                             2.0 for i in range(len(dz1) - 1)])
            unew = np.interp(znew, z, u)
            sm_i.append(unew)
        return np.array(sm_i)

    def fertilizers(self, planting,lon,lat):
        """Setup fertilizer information."""
        if not self.fertilizer or self.fertilizer == "default":
            fertilizers = None
        else:
            fertilizers = self.readFertilizer(planting,lon,lat)
        return fertilizers

    def _ap(self):
        """Return code for fertilizer of specific crop."""
        ap = dict(rice="014", maize="001")
        if self.crop in ap:
            return ap[self.crop]
        else:
            return "014"

    def readFertilizer(self, planting,lon,lat):
        """Extract fertilizer information from the database."""
        log = logging.getLogger(__name__)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        self.fertilizer = 'dynamic'
        if dbio.tableExists(self.dbname, "dssat", "fertilizer"):                
            fertilizers=[]
            qens = 0
            for ens in range(1):
                fertilizer = []
                qens += 1
                sql = "select lag, amount from dssat.fertilizer as c, {0}.agareas as a where crop='{1}' and st_intersects(c.geom,st_geomfromtext('POINT({2} {3})',4326))".format(self.name, self.crop,lon,lat)
                cur.execute(sql)
                results = cur.fetchall()
                results = sorted(list(set(results)))
                for r in results:
                    fe = (planting+timedelta(r[0]), "005", self._ap(), 1.0, r[1]) 
                    fertilizer.append(fe)
                fertilizers.append(fertilizer)
        else:
            fertilizers = None
            log.warning("Default fertilizer used despite user dynamic option.")
        cur.close()
        db.close()
        print(fertilizers)
        #sys.exit()
        return fertilizers

    def Irrigation(self, startdate,lon,lat):
        """Setup fertilizer information."""
        if not self.irrigation or self.irrigation == "default":
            irrigation = None
        else:
            irrigation = self.readIrrigation(startdate,lon,lat)
        return irrigation

    def readIrrigation(self, startdate,lon,lat):
        """Extract irrigation information from the database."""
        log = logging.getLogger(__name__)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        if dbio.tableExists(self.dbname, "dssat", "irrigation"):
            sql = "select lag from dssat.irrigation as c where st_intersects(c.geom,st_geomfromtext('POINT({0} {1})',4326))".format(lon,lat)
            cur.execute(sql)
            if bool(cur.rowcount):
                cur.execute("select max(ensemble) from dssat.irrigation as c where st_intersects(c.geom,st_geomfromtext('POINT({0} {1})',4326))".format(lon,lat))
                max_nens = cur.fetchone()
                cur.execute("select max(fdate) from dssat.irrigation as c where st_intersects(c.geom,st_geomfromtext('POINT({0} {1})',4326))".format(lon,lat))
                max_fdate = cur.fetchone()
                print('max ensemble',max_nens)
                
                Irrigation=[]
                qens = 0
                for ens in range(self.nens):
                    irrigation = []
                    if qens > int(max_nens[0]):
                        qens = 0
                    qens += 1
                    sql = "select lag, amount from dssat.irrigation as c, {0}.agareas as a where ensemble='{1}' and crop='{2}' and fdate='{3}' and st_intersects(c.geom,st_geomfromtext('POINT({4} {5})',4326))".format(self.name, qens, self.crop,max_fdate[0],lon,lat)
                    cur.execute(sql)
                    results = cur.fetchall()
                    ir=0
                    for r in results:
                        IR = (startdate[qens-1]+timedelta(r[0]), "IR0{0}".format(str(ir).rjust(2, '0')), r[1]) 
                        ir+=1
                        irrigation.append(IR)
                    Irrigation.append(irrigation)               
            else:
                Irrigation = None
                log.warning("Default irrigation used because the geom of the area is not foung")
        else:
            Irrigation = None
            log.warning("Default irrigation used despite user dynamic option.")
        cur.close()
        db.close()
        return Irrigation


    def copyModelFiles(self, geom, pi, dssatexe):
        """Copy DSSAT model files to instance's directory."""
        gid, lat, lon = geom
        modelpath = os.path.abspath("{0}/{1}_{2}_{3}".format(self.path, lat, lon, pi))
        self.modelpaths[(gid, pi)] = modelpath
        os.mkdir(modelpath)
        os.mkdir(modelpath + "/ENKF_Results")
        dssatpath = os.path.dirname(os.path.abspath(__file__)) + "/../../external/dssat"
        shutil.copyfile("{0}/{1}".format(dssatpath, dssatexe), "{0}/{1}".format(modelpath, dssatexe))
        # ensure that executable has right permissions
        os.chmod("{0}/{1}".format(modelpath, dssatexe), 0o100)
        distutils.dir_util.copy_tree("{0}/dssat".format(self.data_path), modelpath)

    def setupModelInstance(self, geom, dssatexe):
        """Setup parameters and write input files for a DSSAT model instance
        over a specific geometry."""
        log = logging.getLogger(__name__)
        gid, lon, lat = geom
        c = np.argmin(np.sqrt((lat - self.lat) **
                              2 + (lon - self.lon) ** 2))
        # use the soil depths from the nearest VIC pixel to the centroid
        depths = np.array(self.depths[c])
        planting = self.planting(lat, lon)
        # FIXME: instead of 180 days perhaps use information from harvest dates to define the simulation end date?
        harvest_days = 180
        # we will start the simulation 7 days before the first planting date if possible
        simstartdt = planting[0] - timedelta(7)
        year, month, day, weather, sm, vlai = self.readVICOutput(gid, depths, simstartdt, planting[-1] + timedelta(harvest_days))
        vicstartdt = date(year[0], month[0], day[0])
        for pi, pdt in enumerate(planting):
            if (pdt - vicstartdt).days < 0:
                log.warning("Cannot perform simulation for planting date {0}. Earliest available VIC output is on {1}".format(pdt.strftime("%Y-%m-%d"), vicstartdt.strftime("%Y-%m-%d")))
            else:
                simstartdt = max(vicstartdt, pdt - timedelta(7))
                self.copyModelFiles(geom, pi, dssatexe)
                modelpath = self.modelpaths[(gid, pi)]
                self.modelstart[(gid, pi)] = simstartdt
                fertilizers=self.fertilizers(planting[pi],lon,lat)
                irrigation=self.Irrigation(simstartdt,lon,lat)
                print(fertilizers,irrigation)
                dz, smi = self.writeControlFile(modelpath, sm, depths, simstartdt, gid, self.lat[c], self.lon[c], pdt, fertilizers, irrigation)
                ti0 = [i for i in range(len(year)) if simstartdt == date(year[i], month[i], day[i])][0]
                ti1 = ti0 + harvest_days
                if ti1 > len(year):
                    log.warning("Inadequate record legnth in VIC data to ensure harvest for {0} planting date! Plant will not reach maturity and yield values will be invalid. Please extent VIC simulation to at least {1}!".format(pdt.strftime("%Y-%m-%d"), (pdt + timedelta(harvest_days)).strftime("%Y-%m-%d")))
                    ti1 = len(year) - 1
                self.writeWeatherFiles(modelpath, self.name, year, month, day, weather, self.elev[c], self.lat[c], self.lon[c], ti0, ti1)
                self.writeSoilMoist(modelpath, year, month, day, smi, dz, ti0, ti1)
                self.writeLAI(modelpath, gid, year, month, day, ti0, ti1, viclai=vlai)
                self.writeConfigFile(modelpath, smi.shape[1], simstartdt, date(year[ti1], month[ti1], day[ti1]))
                log.info("Wrote DSSAT for planting date {0}".format(pdt.strftime("%Y-%m-%d")))

    def runModelInstance(self, modelpath, dssatexe):
        """Runs DSSAT model instance."""
        log = logging.getLogger(__name__)
        os.chdir(modelpath)
        crop = {'maize': 'MZ', 'rice': 'RI'}
        if bool(self.assimilate):
            if str(self.assimilate).lower() == "sm":
                sm_assim = "Y"
                lai_assim = "N"
            elif str(self.assimilate).lower() == "lai":
                sm_assim = "N"
                lai_assim = "Y"
            else:
                sm_assim = lai_assim = "Y"
        else:
            sm_assim = lai_assim = "N"
        proc = subprocess.Popen(["./{}".format(dssatexe), "SOIL_MOISTURE.ASC", "LAI.txt", "SM{0}".format(sm_assim), "LAI{0}".format(lai_assim), "{0}".format(crop[self.crop])])
        out, err = proc.communicate()
        log.debug(out)

    def save(self):
	# this section saves the dssat output at each date (includes fdate - old version)        
	# every run it deletes the old schema of the same name is present 
        # does not add new data to already present schema
        """Saves DSSAT output to database."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        #cur.execute("DROP TABLE IF EXISTS {0}.dssat_all".format(self.name))
        #if not bool(cur.rowcount):
        #cur.execute("create table {0}.dssat_all (id serial primary key, gid int, ensemble int, fdate date, wsgd real, lai real, gwad real, geom geometry, CONSTRAINT enforce_dims_geom CHECK (st_ndims(geom) = 2), CONSTRAINT enforce_geotype_geom CHECK (geometrytype(geom) = 'POLYGON'::text OR geometrytype(geom) = 'MULTIPOLYGON'::text OR geom IS NULL))".format(self.name))
        cur.execute("create table if not exists {0}.dssat_all (id serial primary key, gid int, season varchar(10), ensemble int, fdate date, wsgd real, lai real, gwad real)".format(self.name))
        db.commit()
        
	# overwrite overlapping dates
        cur.execute("delete from {0}.dssat_all where fdate>=date'{1}-{2}-{3}' and fdate<=date'{4}-{5}-{6}'".format(self.name, self.startyear, self.startmonth, self.startday, self.endyear, self.endmonth, self.endday))
        sql = "insert into {0}.dssat_all (fdate, gid, season, ensemble, gwad, wsgd, lai) values (%(dt)s, %(gid)s, %(season)s, %(ens)s, %(gwad)s, %(wsgd)s, %(lai)s)".format(self.name)
        for gid, pi in self.modelpaths:
            modelpath = self.modelpaths[(gid, pi)]
            startdt = self.modelstart[(gid, pi)]
            # get season (season 1 or 2) from csv to aggregate yields by season instead of planting
            # NOTE : as more countries are added they must be added to csv, and double check the logic
            seasondf = pd.read_csv('/data/RHEAS/data/crops/calendar/maize_season.csv')
            seasondf = seasondf.loc[seasondf['country']==self.dbname]
            if(len(seasondf)==1):
                season = str(startdt.year)+' 1'
            else:
                jsdt = startdt.timetuple().tm_yday
                pstart = seasondf['planting start'].loc[seasondf['season']==2].values[0]
                pend = seasondf['planting end'].loc[seasondf['season']==2].values[0]
                if((pstart < jsdt)&(jsdt < pend)):
                    season = str(startdt.year)+' 2'
                elif(startdt.month==12):
                    season = str(int(startdt.year) + 1)+' 1'
                else:
                    season = str(startdt.year)+' 1'
            print(startdt,season)
            for e in range(self.nens):
                if os.path.exists("{0}/PLANTGRO{1:03d}.OUT".format(modelpath, e + 1)):
                    with open("{0}/PLANTGRO{1:03d}.OUT".format(modelpath, e + 1)) as fin:
                        line = fin.readline()
                        while line.find("YEAR") < 0:
                            line = fin.readline()
                        for line in fin:
                            data = line.split()
                            dt = date(int(data[0]), 1, 1) + \
                                timedelta(int(data[1]) - 1)
                            dts = "{0}-{1}-{2}".format(dt.year, dt.month, dt.day)
                            if self.cultivars[gid][e] is None:
                                cultivar = ""
                            else:
                                cultivar = self.cultivars[gid][e]
                            #if float(data[9]) > 0.0:
                            cur.execute(sql, {'dt': dts, 'season': season, 'ens': e + 1, 'gwad': float(
                                    data[9]), 'wsgd': float(data[18]), 'lai': float(data[6]), 'gid': gid, 'cultivar': cultivar})
        #cur.execute("update {0}.dssat_all as d set geom = a.geom from {0}.agareas as a where a.gid=d.gid".format(self.name))
        db.commit()
        cur.execute("drop index if exists {0}.d_t".format(self.name))
        #cur.execute("drop index if exists {0}.d_s".format(self.name))
        cur.execute("create index d_t on {0}.dssat_all(fdate)".format(self.name))
        #cur.execute("create index d_s on {0}.dssat_all using gist(geom)".format(self.name))
        db.commit()
        cur.close()
        db.close()
        #self.yieldTable()


	# this section only saves the values at harvest date
	# Modified schema to include county name and code ** attribute names are hard coded **
	

        """Saves DSSAT output to database."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute(
            "select * from information_schema.tables where table_name='dssat' and table_schema='{0}'".format(self.name))
        if not bool(cur.rowcount):
            cur.execute("create table {0}.dssat (id serial primary key, gid int, season varchar(10), cname varchar(50), ccode varchar(20), ensemble int, harvest date, planting date, wsgd real, lai real, gwad real, geom geometry, CONSTRAINT enforce_dims_geom CHECK (st_ndims(geom) = 2), CONSTRAINT enforce_geotype_geom CHECK (geometrytype(geom) = 'POLYGON'::text OR geometrytype(geom) = 'MULTIPOLYGON'::text OR geom IS NULL))".format(self.name))
            #cur.execute("create table {0}.dssat (id serial primary key, gid int, ensemble int, harvest date, planting date, wsgd real, lai real, gwad real, geom geometry, CONSTRAINT enforce_dims_geom CHECK (st_ndims(geom) = 2), CONSTRAINT enforce_geotype_geom CHECK (geometrytype(geom) = 'POLYGON'::text OR geometrytype(geom) = 'MULTIPOLYGON'::text OR geom IS NULL))".format(self.name))
            db.commit()
        # overwrite overlapping dates
        cur.execute("delete from {0}.dssat where planting>=date'{1}-{2}-{3}' and planting<=date'{4}-{5}-{6}'".format(self.name, self.startyear, self.startmonth, self.startday, self.endyear, self.endmonth, self.endday))
        sql = "insert into {0}.dssat (planting, harvest, gid, season, ensemble, gwad, wsgd, lai) values (%(pdt)s, %(hdt)s, %(gid)s, %(season)s, %(ens)s, %(gwad)s, %(wsgd)s, %(lai)s)".format(self.name)
        for gid, pi in self.modelpaths:
            modelpath = self.modelpaths[(gid, pi)]
            startdt = self.modelstart[(gid, pi)]
            # get season (season 1 or 2) from csv to aggregate yields by season instead of planting
            # NOTE : as more countries are added they must be added to csv, and double check the logic
            seasondf = pd.read_csv('/data/RHEAS/data/crops/calendar/maize_season.csv')
            seasondf = seasondf.loc[seasondf['country']==self.dbname]
            if(len(seasondf)==1):
                season = str(startdt.year)+' 1'
            else:
                jsdt = startdt.timetuple().tm_yday
                pstart = seasondf['planting start'].loc[seasondf['season']==2].values[0]
                pend = seasondf['planting end'].loc[seasondf['season']==2].values[0]
                if((pstart < jsdt)&(jsdt < pend)):
                    season = str(startdt.year)+' 2'
                elif(startdt.month==12):
                    season = str(int(startdt.year) + 1)+' 1'
                else:
                    season = str(startdt.year)+' 1'

            for e in range(self.nens):
                try: 
                    data = pd.read_csv("{0}/PLANTGRO{1:03d}.OUT".format(modelpath, e + 1), delim_whitespace=True)
                    # FIXME: Should we write the entire time series here?
                    hvst = data.index[-1]  # assume that harvest occurred at last date of simulation
                    plnt = data.index[0]
                    year = data['@YEAR']
                    doy = data['DOY']
                    lai = data['LAID'][hvst]
                    gwad = data['GWAD'][hvst].astype(float)
                    wsgd = data['WSGD'][hvst]
                    hdts = datetime.strptime("{0:04d}{1:03d}".format(year[hvst], doy[hvst]), "%Y%j").strftime("%Y-%m-%d")
                    pdts = datetime.strptime("{0:04d}{1:03d}".format(year[plnt], doy[plnt]), "%Y%j").strftime("%Y-%m-%d")
                    if self.cultivars[gid][e] is None:
                        cultivar = ""
                    else:
                        cultivar = self.cultivars[gid][e]
                    if gwad > 0:
                        cur.execute(sql, {'pdt': pdts, 'hdt': hdts, 'season': season, 'ens': e + 1, 'gwad': gwad, 'wsgd': wsgd, 'lai': lai, 'gid': gid, 'cultivar': cultivar})
                except:
                    print('Could not save output from file {0}/PLANTGRO{1:03d}.OUT. Continuing...'.format(modelpath, e + 1))
        cur.execute(
            "update {0}.dssat as d set geom = a.geom from {0}.agareas as a where a.gid=d.gid".format(self.name))
        db.commit()
        cur.execute("update  {0}.dssat as d set cname = a.NAME1 from {0}.agareas as a where a.gid = d.gid".format(self.name))
        cur.execute("update  {0}.dssat as d set ccode = a.GID1 from {0}.agareas as a where a.gid = d.gid".format(self.name))
        db.commit()
        cur.execute("drop index if exists {0}.d_t".format(self.name))
        cur.execute("drop index if exists {0}.d_s".format(self.name))
        cur.execute(
            "create index d_t on {0}.dssat(planting)".format(self.name))
        cur.execute(
            "create index d_s on {0}.dssat using gist(geom)".format(self.name))
        db.commit()
        cur.close()
        db.close()
        self.yieldTable()
        #change folder name to dbname+schemaname
        changename = 'mv /data/RHEAS/{0} /data/RHEAS/{1}{2}'.format(self.path,self.dbname,self.name)
        os.system(changename)


    def yieldTable(self):
        """Create table for crop yield statistics."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        if dbio.tableExists(self.dbname, self.name, "yield"):
            cur.execute("drop table {0}.yield".format(self.name))
        sql = "create table {0}.yield as (select gid, geom, cname, ccode, season, min(planting) as planting, min(harvest) as first_harvest, max(harvest) as last_harvest, avg(gwad) as avg_yield, max(gwad) as max_yield, min(gwad) as min_yield, stddev(gwad) as std_yield from {0}.dssat group by gid,geom,cname,ccode,season)".format(self.name)
        cur.execute(sql)
        cur.execute("alter table {0}.yield add column crop text".format(self.name))
        # add production column for zambia based on csv of cropland area per county
        if (self.dbname == 'zambia'):
            area = pd.read_csv('/data/RHEAS/data/crops/mask/zambia_maize.csv')
            cur.execute("alter table {0}.yield add column production real".format(self.name))
            for c in area.NAME1.unique():
                #skip if county doesn't have any yield in db
                c_area = area['Area_ha'].loc[area['NAME1']==c].iloc[0]
                sql = "SELECT exists (SELECT cname FROM zmb_n_25.yield WHERE cname = \'{0}\' LIMIT 1)".format(c)
                cur.execute(sql)
                result = cur.fetchone()
                if (result[0] == True):
                    cur.execute("update {0}.yield set production = (avg_yield * {1}) where cname = \'{2}\'".format(self.name,c_area,c))
        cur.execute("drop index if exists {0}.yield_s".format(self.name))
        db.commit()
        cur.execute("create index yield_s on {0}.yield using gist(geom)".format(self.name))
        cur.close()
        db.close()

    def run(self, dssatexe="mdssat", crop_threshold=0.05):
        """Runs DSSAT simulation."""
        self.readVICSoil()
        geoms = self.readShapefile()
        cropfract = self.calcCroplandFract()
        for geom in geoms:
            gid = geom[0]
            if cropfract[gid] >= crop_threshold:
                self.setupModelInstance(geom, dssatexe)
        pwd = os.getcwd()
        for k in self.modelpaths:
            modelpath = self.modelpaths[k]
            self.runModelInstance(modelpath, dssatexe)
        self.save()
        os.chdir(pwd)
