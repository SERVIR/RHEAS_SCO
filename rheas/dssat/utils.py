""" Module for helper function of the DSSAT model

.. module:: utils
   :synopsis: Definition of the DSSAT model utility functions

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import os
import random
import string
import subprocess

from .. import dbio


def addCultivar(dbname, shapefile, params, nens=40, crop="maize"):
    """Add cultivar parameters to the database *dbname* corresponding
    to the area defined in the *shapefile*. The *params* is a list of dictionaries,
    where the keys of each dictionary correspond to parameters, and each object in
    the list corresponds to a cultivar variant. The *nens* parameters is the size
    of the ensemble to be created."""
    temptable = ''.join(random.SystemRandom().choice(
        string.ascii_letters) for _ in range(8))
    if os.path.exists(shapefile):
        subprocess.call("shp2pgsql -d -s 4326 -g geom {0} {1} | psql -d {2}".format(
            shapefile, temptable, dbname), shell=True)
        db = dbio.connect(dbname)
        cur = db.cursor()
        e = 0
        while e < nens:
            for c in range(len(params)):
                if crop == "maize" and all(p in params[c] for p in ['p1', 'p2', 'p5', 'g2', 'g3', 'phint']):
                    if e < nens:
                        sql = "insert into dssat.cultivars (geom) (select geom from {0})".format(
                            temptable)
                        cur.execute(sql)
                        sql = "update dssat.cultivars set crop='maize',ensemble={0},{1} where ensemble is null".format(
                            e + 1, ",".join(["{0}={1}".format(k, params[c][k]) for k in params[c]]))
                        cur.execute(sql)
                        e += 1
                elif crop == "rice" and all(p in params[c] for p in ['p1', 'p2r', 'p5', 'p2o', 'g1', 'g2', 'g3', 'g4']):
                    if e < nens:
                        sql = "insert into dssat.cultivars (geom) (select geom from {0})".format(
                            temptable)
                        cur.execute(sql)
                        sql = "update dssat.cultivars set crop='rice',ensemble={0},{1} where ensemble is null".format(
                            e + 1, ",".join(["{0}={1}".format(k, params[c][k]) for k in params[c]]))
                        cur.execute(sql)
                        e += 1
                else:
                    print("Missing parameters for {0} crop".format(crop))
                    params.pop(c)  # remove element with missing parameters
                    break
        cur.execute("drop table {0}".format(temptable))
        db.commit()
        cur.close()
        db.close()
    else:
        print(
            "Shapefile {0} cannot be found. Not adding cultivars!".format(shapefile))



def addfertilizer(dbname, shapefile, params, nens,dt,name,crop):
    """Add fertilizer amount to the database *dbname* corresponding
    to the area defined in the *shapefile*. The *params* is a list of dictionaries,
    where the keys of each dictionary correspond to the number of days after planting, and each object in
    the list corresponds to the amount of fertilizer (kg/ha). The *nens* parameters is the size
    of the ensemble to be created."""
    DOY=list(params.keys())

    FER=list(params.values())

    

    temptable = ''.join(random.SystemRandom().choice(string.ascii_letters) for _ in range(8))
    #temptable='dssat.'+a
    
    if os.path.exists(shapefile):
        subprocess.call("shp2pgsql -d -s 4326 -g geom {0} {1} | psql -d {2}".format(
            shapefile, temptable, dbname), shell=True)

        db = dbio.connect(dbname)
        cur = db.cursor()
        cur.execute("alter table {0} add column fdate text, add column ensemble text, add column name text, add column crop text, add column lag float, add column amount float".format(temptable))
        if dbio.tableExists(dbname, "dssat", "fertilizer"):    
            for m in range(len(params)):
                print(m)
                e=0
                while e < nens:
                    print(e)
                    cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', name='{5}', crop='{6}'".format(temptable,e+1,DOY[m],FER[m],dt,name,crop))                       
                    cur.execute("insert into dssat.fertilizer (gid,fdate,geom,ensemble,name,crop,lag,amount) select gid,fdate,geom,ensemble,name,crop,lag,amount from {0}".format(temptable))
                    e += 1
        else:
            cur.execute("create table dssat.fertilizer (gid int,fdate text, geom geometry,ensemble text, name text, crop text, lag float, amount float);")         
            
            for m in range(len(params)):
                e = 0
                while e < nens:
                    if e < nens:
                        
                        cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', name='{5}', crop='{6}'".format(temptable,e+1,DOY[m],FER[m],dt,name,crop))                       
                        cur.execute("insert into dssat.fertilizer (gid,fdate,geom,ensemble,name,crop,lag,amount) select gid,fdate,geom,ensemble,name,crop,lag,amount from {0}".format(temptable))

                        e += 1
                    else:
                        print("Missing parameters")
                        # remove element with missing parameters
                        break
        
        cur.execute("truncate table {0}".format(temptable))
        db.commit()
        cur.close()
        db.close()
    else:
        print("Shapefile {0} cannot be found. Not adding cultivars!".format(shapefile))




def addirrigation(dbname, shapefile, params, nens,dt, IR_pr,crop):
    """Add irrigation amount to the database *dbname* corresponding
    to the area defined in the *shapefile*. The *params* is a list of dictionaries,
    where the keys of each dictionary correspond to the number of days after planting, and each object in
    the list corresponds to the amount of irrigation (mm). The *nens* parameters is the size
    of the ensemble to be created."""
    DOY=list(params.keys())

    IR=list(params.values())   

    temptable = ''.join(random.SystemRandom().choice(string.ascii_letters) for _ in range(8))
    #temptable='dssat.'+a
    
    if os.path.exists(shapefile):
        subprocess.call("shp2pgsql -d -s 4326 -g geom {0} {1} | psql -d {2}".format(
            shapefile, temptable, dbname), shell=True)
        db = dbio.connect(dbname)
        cur = db.cursor()
        cur.execute("alter table {0} add column fdate text, add column ensemble text, add column crop text, add column lag float, add column amount float".format(temptable))
        if dbio.tableExists(dbname, "dssat", "irrigation"):    
            for m in range(len(params)):
                print(m)
                e=0
                while e <int(nens*IR_pr/100):
                    #print(e)
                    cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', crop='{5}'".format(temptable,e+1,DOY[m],IR[m],dt,crop))                       
                    cur.execute("insert into dssat.irrigation (gid,fdate,geom,ensemble,crop,lag,amount) select gid,fdate,geom,ensemble,crop,lag,amount from {0}".format(temptable))
                    e += 1
                e=0
                while e <nens - int(nens*IR_pr/100):
                    #print(e)
                    cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', crop='{5}'".format(temptable,e+1+int(nens*IR_pr/100),0.0,0.0,dt,crop))                       
                    cur.execute("insert into dssat.irrigation (gid,fdate,geom,ensemble,crop,lag,amount) select gid,fdate,geom,ensemble,crop,lag,amount from {0}".format(temptable))
                    e +=1
        else:
            cur.execute("create table dssat.irrigation (gid int,fdate text, geom geometry,ensemble text, crop text, lag float, amount float);")         
            for m in range(len(params)):
                e=0
                while e < int(nens*IR_pr/100):
                    #print(e)
                    cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', crop='{5}'".format(temptable,e+1,DOY[m],IR[m],dt,crop))                       
                    cur.execute("insert into dssat.irrigation (gid,fdate,geom,ensemble,crop,lag,amount) select gid,fdate,geom,ensemble,crop,lag,amount from {0}".format(temptable))
                    e += 1
                e=0
                while e<nens-int(nens*IR_pr/100):
                    #print(e)
                    cur.execute("update {0} set ensemble = '{1}',lag = '{2}',amount='{3}',fdate='{4}', crop='{5}'".format(temptable,e+1+int(nens*IR_pr/100),0.0,0.0,dt,crop))                       
                    cur.execute("insert into dssat.irrigation (gid,fdate,geom,ensemble,crop,lag,amount) select gid,fdate,geom,ensemble,crop,lag,amount from {0}".format(temptable))
                    e +=1
                    
        
        cur.execute("truncate table {0}".format(temptable))
        db.commit()
        cur.close()
        db.close()
    else:
        print("Shapefile {0} cannot be found. Not adding cultivars!".format(shapefile))
