#!/usr/bin/python
import logging
import os, sys
import argparse
import shutil
import zipfile
import math
import numpy as np
import osr
import glob
import saa_func_lib as saa
from byteSigmaScale import byteSigmaScale
from createAmp import createAmp
from get_zone import get_zone
from getSubSwath import get_bounding_box_file
from execute import execute
from getParameter import getParameter
from makeAsfBrowse import makeAsfBrowse
from osgeo import gdal

# Convert point from geographic to UTM projection
def geometry_geo2proj(lat_max,lat_min,lon_max,lon_min):
    zone = get_zone(lon_min,lon_max) 
    if (lat_min+lat_max)/2 > 0:
        proj = ('326%02d' % int(zone))
    else:
        proj = ('327%02d' % int(zone))
	
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(4326)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(int(proj))
    coordTrans = osr.CoordinateTransformation(inSpatialRef,outSpatialRef)
  
    x1, y1, h = coordTrans.TransformPoint(lon_max, lat_min)
    logging.debug("Output coordinate: {} {} {}".format(x1,y1,h))
    x2, y2, h = coordTrans.TransformPoint(lon_min, lat_min)
    logging.debug("Output coordinate: {} {} {}".format(x2,y2,h))
    x3, y3, h = coordTrans.TransformPoint(lon_max, lat_max)
    logging.debug("Output coordinate: {} {} {}".format(x3,y3,h))
    x4, y4, h = coordTrans.TransformPoint(lon_min, lat_max)
    logging.debug("Output coordinate: {} {} {}".format(x4,y4,h))

    y_min = min(y1,y2,y3,y4)
    y_max = max(y1,y2,y3,y4)
    x_min = min(x1,x2,x3,x4)
    x_max = max(x1,x2,x3,x4)
    
    false_easting = outSpatialRef.GetProjParm(osr.SRS_PP_FALSE_EASTING)
    false_northing = outSpatialRef.GetProjParm(osr.SRS_PP_FALSE_NORTHING)

    return zone, false_northing, y_min, y_max, x_min, x_max


def create_dem_par(basename,dataType,pixel_size,lat_max,lat_min,lon_max,lon_min):
    demParIn = "dem_par.in"
    zone, false_north, y_min, y_max, x_min, x_max = geometry_geo2proj(lat_max,lat_min,lon_max,lon_min)

    logging.debug("Output Coordinates: {} {} {} {}".format(y_min, y_max, x_min, x_max))
    
    f = open(demParIn,"w")
    f.write("UTM\n")
    f.write("WGS84\n")
    f.write("1\n")
    f.write("{}\n".format(zone))
    f.write("{}\n".format(false_north))
    f.write("{}\n".format(basename))
    if "float" in dataType:
        f.write("REAL*4\n")
    elif "int16" in dataType:
        f.write("INTEGER*2\n")
    f.write("0.0\n")
    f.write("1.0\n")
    
    xsize = np.floor(abs((x_max-x_min)/pixel_size))
    ysize = np.floor(abs((y_max-y_min)/pixel_size))
    
    f.write("{}\n".format(int(xsize)))
    f.write("{}\n".format(int(ysize)))
    f.write("{} {}\n".format(-1.0*pixel_size,pixel_size))
    f.write("{} {}\n".format(y_max,x_min))
    f.close()
    
    return demParIn

def blank_bad_data(rawFile,x,y,left=15,right=15):

    # Read in the data
    data = np.fromfile(rawFile,dtype=np.float32)
    data = np.reshape(data,(y,x))
    data = data.byteswap()

    # For each line in the file
    for i in xrange(y):
        # Black out the start of the line 
        for j in xrange(x):
            if data[i,j] != 0:
                data[i,:j+left] = 0
                break
        # Black out the end of the line 
        for j in xrange(x-1,0,-1):
            if data[i,j] != 0:
                data[i,j-right:] = 0
                break

    # Write out the data            
    data = data.byteswap()
    data.tofile(rawFile)

def getFileType(myfile):
    if "GRD" in myfile:
        type = "GRD"
    else:
        type = "SLC"
    return(type)



def geocode_sentinel(infile,outfile,pixel_size=30.0,height=0):

    if not os.path.exists(infile):
        logging.error("ERROR: Input file {} does not exist".format(infile))
        exit(1)
    if "zip" in infile:
        zip_ref = zipfile.ZipFile(myfile, 'r')
        zip_ref.extractall(".")
        zip_ref.close()    
        infile = infile.replace(".zip",".SAFE")

    type = getFileType(infile)

    # Create par file covering the area we want to geocode
    lat_max,lat_min,lon_max,lon_min = get_bounding_box_file(infile)
    logging.debug("Input Coordinates: {} {} {} {}".format(lat_max,lat_min,lon_max,lon_min))
    demParIn = create_dem_par("area_map","float",pixel_size,lat_max,lat_min,lon_max,lon_min)
    execute("create_dem_par area_map.par < {}".format(demParIn),uselogging=True)
    
    # Try to get the precision state vectors
    try:
        cmd = "get_orb.py {}".format(infile)
        logging.debug("Getting precision orbit information")
        execute(cmd,uselogging=True)
    except:
        logging.warning("Unable to fetch precision state vectors... continuing")
    
    vvcnt = 0
    vvlist = glob.glob("{}/*/*vv*.tiff".format(infile))
    vhlist = glob.glob("{}/*/*vh*.tiff".format(infile))
    hhlist = glob.glob("{}/*/*hh*.tiff".format(infile))
    hvlist = glob.glob("{}/*/*hv*.tiff".format(infile))
    
    print vvlist, vhlist, hhlist, hvlist
    
    if vvlist:
        pol = "vv"
        process_pol(pol,type,infile,outfile,pixel_size,height)
        if vhlist:
             process_pol("vh",type,infile,outfile,pixel_size,height)     
    if hhlist:
        pol = "hh"
        process_pol(pol,type,infile,outfile,pixel_size,height) 
        if hvlist:
            process_pol("hv",type,infile,outfile,pixel_size,height)   
	
    # Create geotiff and ASF browse images
    basename = "{out}_{pol}".format(out=outfile,pol=pol)
    tiffile = "{}.tif".format(basename)

    ampfile = createAmp(tiffile,nodata=0)
    newfile = ampfile.replace(".tif","_sigma.tif")
    byteSigmaScale(ampfile,newfile)
    os.remove(ampfile)
    makeAsfBrowse(newfile,outfile)
    os.remove(newfile)
    os.mkdir("PRODUCT")
    shutil.move(tiffile,"PRODUCT")
    shutil.move("{}.png.aux.xml".format(outfile),"PRODUCT")
    shutil.move("{}_large.png.aux.xml".format(outfile),"PRODUCT")
    shutil.move("{}.png".format(outfile),"PRODUCT")
    shutil.move("{}_large.png".format(outfile),"PRODUCT")
    shutil.move("{}.kmz".format(outfile),"PRODUCT")
    


def process_pol(pol,type,infile,outfile,pixel_size,height):

    logging.info("Processing the {} polarization".format(pol))

    grd = "{out}.{pol}.grd".format(out=outfile,pol=pol)
    mgrd = "{out}.{pol}.mgrd".format(out=outfile,pol=pol)
    utm = "{out}.{pol}.utm".format(out=outfile,pol=pol)

    # Ingest the granule into gamma format
    if "GRD" in type:
        cmd = "par_S1_GRD {inf}/*/*{pol}*.tiff {inf}/*/*{pol}*.xml {inf}/*/*/calibration-*{pol}*.xml \
              {inf}/*/*/noise-*{pol}*.xml {grd}.par {grd}".format(inf=infile,pol=pol,grd=grd)
        execute(cmd,uselogging=True)
	
        # Update the state vectors
        try:
            for eoffile in glob.glob("*.EOF"):
                logging.debug("Applying precision orbit information")
                cmd = "S1_OPOD_vec {grd}.par {eof}".format(grd=grd,eof=eoffile)
                execute(cmd,uselogging=True)
        except:
            logging.warning("Unable to get precision state vectors... continuing...")

        # Multi-look the image
        look_fact = np.floor((pixel_size/10.0)+0.5) 
        if look_fact > 1.0:
            cmd = "multi_look_MLI {grd} {grd}.par {mgrd} {mgrd}.par {lks} {lks}".format(grd=grd,mgrd=mgrd,lks=look_fact)
            execute(cmd,uselogging=True)
        else:
	    shutil.copy(grd,mgrd)
            shutil.copy("{}.par".format(grd),"{}.par".format(mgrd))

    else:
        par_s1_slc(pol)
	burst_tab = get_bursts(infile)
  

    # Blank out the bad data at the left and right edges
    dsx = int(getParameter("{}.par".format(mgrd),"range_samples",uselogging=True))
    dsy = int(getParameter("{}.par".format(mgrd),"azimuth_lines",uselogging=True))
    blank_bad_data(mgrd, dsx, dsy, left=40, right=8)


    # Create geocoding look up table
    # by running either gec_map or gec_map_grd

    #
    # This command doesn't work, telling me no DEM overlap in longitude/easting!
#    cmd = "gec_map_grd {mgrd}.par area_map.par {height} small_map.par small_map.utm_to_rdc - -".format(mgrd=mgrd,height=height)
#    execute(cmd,uselogging=True)

    cmd = "gec_map {mgrd}.par - area_map.par {height} small_map.par small_map.utm_to_rdc".format(mgrd=mgrd,height=height)
    execute(cmd,uselogging=True)
   
    # Gecode the granule by running geocode_back
    outSize = getParameter("small_map.par","width",uselogging=True)
    cmd = "geocode_back {mgrd} {dsx} small_map.utm_to_rdc {utm} {os}".format(mgrd=mgrd,utm=utm,dsx=dsx,os=outSize)
    execute(cmd,uselogging=True)

    # Create the geotiff file
    tiffile = "{out}_{pol}.tif".format(out=outfile,pol=pol)
    cmd = "data2geotiff small_map.par {utm} 2 {tif}".format(utm=utm,tif=tiffile)
    execute(cmd,uselogging=True)


###########################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='geocode_sentinel',
    description='Geocode a sentinel-1 granule usig Gamma software')
  parser.add_argument("infile",help="Input zip file or SAFE directory")
  parser.add_argument("outfile",help="Name of output geocoded file")
  parser.add_argument("-t","--terrain_height",help="Average terrain height for geocoding",type=float,default=0.0)
  parser.add_argument("-p","--pixel_size",help="Average terrain height for geocoding",type=float,default=30.0)
  args = parser.parse_args()

  logFile = "geocode_sentinel_{}_log.txt".format(os.getpid())
  logging.basicConfig(filename=logFile,format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG)
  logging.getLogger().addHandler(logging.StreamHandler())
  logging.info("Starting run")

  geocode_sentinel(args.infile,args.outfile,height=args.terrain_height,pixel_size=args.pixel_size)

