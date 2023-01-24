################################
# CODE TO EXTRACT RASTER FROM 
# POSTGRES DATABASE AND SAVE
# AS GEOTIFF
# ---------------
# VIKALP MISHRA
# SERVIR SCO
##############################


import sys, os
import psycopg2
from datetime import timedelta
import datetime as dt
import numpy as np

print(sys.argv)
if len(sys.argv)!=7:
    print('*** Not correct number of arguments supplied ***')
    print('')
    print('Correct Usage:')
    print('python pg2tiff.py <db> <schema> <table> <startdate> <enddate> <outdir>')
    print('')
    print('Exiting now....')
    print('****')
else:

    db = sys.argv[1]                    # database name
    schema = sys.argv[2]                # schema name
    table = sys.argv[3]                 # table name
    hst = 'localhost'                  # host is local if run locally
    usr = 'rheas'                       # user name
    pwd = 'servir@2002'                # password saved for now
    print(db,table)
    startdate = sys.argv[4]
    enddate = sys.argv[5]
    outd = sys.argv[6]

    # setup connection to the database using username and password
    conn = psycopg2.connect("dbname={0} user={1} password={2}".format(db,usr,pwd))
    cur = conn.cursor()

    s = np.array(startdate.split('-')).astype(int)
    e = np.array(enddate.split('-')).astype(int)
    sdate = dt.date(s[0],s[1],s[2])
    edate = dt.date(e[0],e[1],e[2])
    tdelta =dt.timedelta(days=1)
    cdate = sdate
    print(sdate,edate)
    if(table=='soil_moist' or table=='soil_temp'):
        sql = "select rid,layer,fdate from {0}.{1} where fdate>='{2}' and fdate<='{3}'".format(schema,table,sdate,edate)
        cur.execute(sql)
        data = cur.fetchall()
    
        dts = []
        ids = []
        lyrs = []
        for d in data:
            dts.append(d[2])
            lyrs.append(d[1])
            ids.append(d[0])
        print(ids)
        for i in range(len(ids)):
                dtt = dts[i]
                sid = ids[i]
                lyr = lyrs[i]
                print(dtt)
                # using gdal to translate extracted raster to GTiff
                t1 = "gdal_translate -of GTiff "
                #t2 = """ "PG:dbname='{1}' user='{2}' password='{3}' schema='{4}' table='{5}' where='fdate=to_date(\'{6}\','YYYY-MM-DD')' mode=1" """.format(hst,db,usr,pwd,schema,table,dtt)
                
                t2 =f"gdal_translate -of GTiff \"PG:dbname='{db}' user='{usr}' password='{pwd}' schema='{schema}' table='{table}' where='fdate=to_DATE(\\'{dtt}\\',\\'YYYY-MM-DD\\') and layer={lyr}' mode='2'\" {outd}/{schema}_{table}_{dtt}_{lyr}.tiff"
                t3 = outd+'/{0}_{1}_{2}_{3}.tiff'.format(schema,table,dtt,lyr)
                print(t2)
                
                os.system(t2)
                #sys.exit()
                cdate = cdate+tdelta
    else:
        sql = "select rid,fdate from {0}.{1} where fdate>='{2}' and fdate<='{3}'".format(schema,table,sdate,edate)
        cur.execute(sql)
        data = cur.fetchall()
    
        dts = []
        ids = []
        for d in data:
            dts.append(d[1])
            ids.append(d[0])
        print(ids)
        for i in range(len(ids)):
                dtt = dts[i]
                sid = ids[i]
                print(dtt)
                # using gdal to translate extracted raster to GTiff
                t1 = "gdal_translate -of GTiff "
                #t2 = """ "PG:dbname='{1}' user='{2}' password='{3}' schema='{4}' table='{5}' where='fdate=to_date(\'{6}\','YYYY-MM-DD')' mode=1" """.format(hst,db,usr,pwd,schema,table,dtt)
                
                t2 =f"gdal_translate -of GTiff \"PG:dbname='{db}' user='{usr}' password='{pwd}' schema='{schema}' table='{table}' where='fdate=to_DATE(\\'{dtt}\\',\\'YYYY-MM-DD\\')' mode='2'\" {outd}/{schema}_{table}_{dtt}.tiff"
                t3 = outd+'/{0}_{1}_{2}.tiff'.format(schema,table,dtt)
                print(t2)
                
                os.system(t2)
                #sys.exit()
                cdate = cdate+tdelta            

