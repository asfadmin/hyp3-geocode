#!/usr/bin/python
import logging
import os, sys
import argparse
import shutil
import zipfile
import math
import osr
import glob
import datetime
import numpy as np
from osgeo import gdal
from lxml import etree

import saa_func_lib as saa
from byteSigmaScale import byteSigmaScale
from createAmp import createAmp
from get_zone import get_zone
from getSubSwath import get_bounding_box_file
from execute import execute
from getParameter import getParameter
from makeAsfBrowse import makeAsfBrowse
from par_s1_slc_single import par_s1_slc_single
from SLC_copy_S1_fullSW import SLC_copy_S1_fullSW
from rtc2color import rtc2color
from asf_geometry import geometry_geo2proj
from getBursts import getBursts
from make_arc_thumb import pngtothumb

def create_dem_par(basename,dataType,pixel_size,lat_max,lat_min,lon_max,lon_min):
    demParIn = "{}_dem_par.in".format(basename)
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

def process_pol(pol,type,infile,outfile,pixel_size,height,make_tab_flag=True,gamma0_flag=False):

    logging.info("Processing the {} polarization".format(pol))

    grd = "{out}.{pol}.grd".format(out=outfile,pol=pol)
    mgrd = "{out}.{pol}.mgrd".format(out=outfile,pol=pol)
    utm = "{out}.{pol}.utm".format(out=outfile,pol=pol)
    area_map = "{out}_area_map.par".format(out=outfile)
    small_map = "{out}_small_map".format(out=outfile)

    look_fact = np.floor((pixel_size/10.0)+0.5)
    if look_fact < 1:
        look_fact = 1 

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
        if look_fact > 1.0:
            cmd = "multi_look_MLI {grd} {grd}.par {mgrd} {mgrd}.par {lks} {lks}".format(grd=grd,mgrd=mgrd,lks=look_fact)
            execute(cmd,uselogging=True)
        else:
	    shutil.copy(grd,mgrd)
            shutil.copy("{}.par".format(grd),"{}.par".format(mgrd))

    else:
        #  Ingest SLC data files into gamma format
        par_s1_slc_single(infile,pol)
        date = infile[17:25]
        burst_tab = getBursts(infile,make_tab_flag)
        shutil.copy(burst_tab,date)

        # Mosaic the swaths together and copy SLCs over        
        back = os.getcwd()
        os.chdir(date) 
        path = "../"
        rlooks = look_fact*5
        alooks = look_fact 
        SLC_copy_S1_fullSW(path,date,"SLC_TAB",burst_tab,mode=2,raml=rlooks,azml=alooks)
        os.chdir(back)
 
        # Rename files
        name = "{}.mli".format(date)
        shutil.move(name,mgrd)
        name = "{}.mli.par".format(date)
        shutil.move(name,"{}.par".format(mgrd))

    if gamma0_flag:
        # Convert sigma-0 to gamma-0
        cmd = "radcal_MLI {mgrd} {mgrd}.par - {mgrd}.sigma - 0 0 -1".format(mgrd=mgrd)
        execute(cmd,uselogging=True)
        cmd = "radcal_MLI {mgrd}.sigma {mgrd}.par - {mgrd}.gamma - 0 0 2".format(mgrd=mgrd)
        execute(cmd,uselogging=True)
        shutil.move("{}.gamma".format(mgrd),mgrd)

    # Blank out the bad data at the left and right edges
    dsx = int(getParameter("{}.par".format(mgrd),"range_samples",uselogging=True))
    dsy = int(getParameter("{}.par".format(mgrd),"azimuth_lines",uselogging=True))

    if "GRD" in type:
        blank_bad_data(mgrd, dsx, dsy, left=20, right=20)

    # Create geocoding look up table
    cmd = "gec_map {mgrd}.par - {amap} {height} {smap}.par {smap}.utm_to_rdc".format(amap=area_map,mgrd=mgrd,height=height,smap=small_map)
    execute(cmd,uselogging=True)
   
    # Gecode the granule
    outSize = getParameter("{}.par".format(small_map),"width",uselogging=True)
    cmd = "geocode_back {mgrd} {dsx} {smap}.utm_to_rdc {utm} {os}".format(mgrd=mgrd,utm=utm,smap=small_map,dsx=dsx,os=outSize)
    execute(cmd,uselogging=True)

    # Create the geotiff file
    tiffile = "{out}_{pol}.tif".format(out=outfile,pol=pol)
    cmd = "data2geotiff {smap}.par {utm} 2 {tif}".format(smap=small_map,utm=utm,tif=tiffile)
    execute(cmd,uselogging=True)

def create_xml_files(infile,outfile,height,type,gamma0_flag):
    # Create XML metadata files
    cfgdir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "config"))
    back = os.getcwd()
    os.chdir("PRODUCT")
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    time = now.strftime("%H%M%S")
    dt = now.strftime("%Y-%m-%dT%H:%M:%S")
    year = now.year
    encoded_jpg = pngtothumb("{}.png".format(outfile))

    if gamma0_flag:
        power_type = "gamma"
    else:
        power_type = "sigma"
    
    for myfile in glob.glob("*.tif"):
        f = open("{}/GeocodingTemplate.xml".format(cfgdir),"r")
        g = open("{}.xml".format(myfile),"w")
	
	if "vv" in myfile:
	    pol = "vv"
	elif "vh" in myfile:
	    pol = "vh"
	elif "hh" in myfile:
	    pol = "hh"
	elif "hv" in myfile:
	    pol = "hv"
	 
        for line in f:
            line = line.replace("[DATE]",date)
            line = line.replace("[TIME]","{}00".format(time))
            line = line.replace("[DATETIME]",dt)
            line = line.replace("[HEIGHT]","{}".format(height))
            line = line.replace("[YEARPROCESSED]","{}".format(year))
            line = line.replace("[YEARACQUIRED]",infile[17:21])
            line = line.replace("[TYPE]",type)
            line = line.replace("[THUMBNAIL_BINARY_STRING]",encoded_jpg)
            line = line.replace("[POL]",pol)
            line = line.replace("[POWERTYPE]",power_type)
            g.write("{}\n".format(line))
        f.close()
        g.close()

    for myfile in glob.glob("*.png"):
    
        if "rgb" in myfile:
            f = open("{}/GeocodingTemplate_color_png.xml".format(cfgdir),"r")
	    encoded_jpg = pngtothumb("{}_rgb.png".format(outfile))
        else:
            f = open("{}/GeocodingTemplate_grayscale_png.xml".format(cfgdir),"r")
	    encoded_jpg = pngtothumb("{}.png".format(outfile))

        if "large" in myfile:  
            res = "medium"
        else:
            res = "low"

        g = open("{}.xml".format(myfile),"w")
        for line in f:
            line = line.replace("[DATE]",date)
            line = line.replace("[TIME]","{}00".format(time))
            line = line.replace("[DATETIME]",dt)
            line = line.replace("[YEARPROCESSED]","{}".format(year))
            line = line.replace("[YEARACQUIRED]",infile[17:21])
            line = line.replace("[TYPE]",type)
            line = line.replace("[THUMBNAIL_BINARY_STRING]",encoded_jpg)
            line = line.replace("[RES]",res)
            g.write("{}\n".format(line))
        f.close()
        g.close()

    os.chdir(back)

def make_products(outfile,pol,cp=None):

    # Create greyscale geotiff and ASF browse images
    tiffile = "{out}_{pol}.tif".format(out=outfile,pol=pol)
    ampfile = createAmp(tiffile,nodata=0)
    newfile = ampfile.replace(".tif","_sigma.tif")
    byteSigmaScale(ampfile,newfile)
    makeAsfBrowse(newfile,outfile)
    os.remove(newfile)

    # Create color ASF browse images
    if cp is not None:
        if pol == "vv":
            basename = "{}_vh".format(outfile)
        else:
            basename = "{}_hv".format(outfile)
        tiffile2  = "{}.tif".format(basename)
        ampfile2 = createAmp(tiffile2,nodata=0)
        outfile2 = ampfile2.replace(".tif","_rgb.tif")
        threshold = -24 
        rtc2color(ampfile,ampfile2, threshold, outfile2, amp=True, cleanup=True)
        colorname = "{}_rgb".format(outfile)
        makeAsfBrowse(outfile2,colorname)
        os.remove(ampfile2)
        os.remove(outfile2)

    os.remove(ampfile)

    # Move results to the PRODUCT directory
    if not os.path.isdir("PRODUCT"):
        os.mkdir("PRODUCT")
    for tiffile in glob.glob("*.tif"):
        shutil.move(tiffile,"PRODUCT")
    for txtfile in glob.glob("*_log.txt"):
        shutil.move(txtfile,"PRODUCT")
    for pngfile in glob.glob("*.png*"):
        shutil.move(pngfile,"PRODUCT")
    for kmzfile in glob.glob("*.kmz"):
        shutil.move(kmzfile,"PRODUCT")


def geocode_sentinel(infile,outfile,pixel_size=30.0,height=0,gamma0_flag=False):

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
    area_map = "{}_area_map".format(outfile)
    demParIn = create_dem_par(area_map,"float",pixel_size,lat_max,lat_min,lon_max,lon_min)
    execute("create_dem_par {}.par < {}".format(area_map,demParIn),uselogging=True)
    
    # Try to get the precision state vectors
    try:
        cmd = "get_orb.py {}".format(infile)
        logging.info("Getting precision orbit information")
        execute(cmd,uselogging=True)
    except:
        logging.warning("Unable to fetch precision state vectors... continuing")
    
    # Get list of files to process
    vvlist = glob.glob("{}/*/*vv*.tiff".format(infile))
    vhlist = glob.glob("{}/*/*vh*.tiff".format(infile))
    hhlist = glob.glob("{}/*/*hh*.tiff".format(infile))
    hvlist = glob.glob("{}/*/*hv*.tiff".format(infile))
    
    crossPol = None 
    if vvlist:
        pol = "vv"
        process_pol(pol,type,infile,outfile,pixel_size,height,make_tab_flag=True,gamma0_flag=gamma0_flag)
        if vhlist:
            process_pol("vh",type,infile,outfile,pixel_size,height,make_tab_flag=False,gamma0_flag=gamma0_flag)     
            crossPol = "vh"
    if hhlist:
        pol = "hh"
        process_pol(pol,type,infile,outfile,pixel_size,height,make_tab_flag=True,gamma0_flag=gamma0_flag) 
        if hvlist:
            process_pol("hv",type,infile,outfile,pixel_size,height,make_tab_flag=False,gamma0_flag=gamma0_flag)   
            crossPol = "hv"
	
    make_products(outfile,pol,cp=crossPol)
    create_xml_files(infile,outfile,height,type,gamma0_flag)
    


###########################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='geocode_sentinel',
    description='Geocode a sentinel-1 granule usig Gamma software')
  parser.add_argument("infile",help="Input zip file or SAFE directory")
  parser.add_argument("outfile",help="Name of output geocoded file")
  parser.add_argument("-t","--terrain_height",help="Average terrain height for geocoding",type=float,default=0.0)
  parser.add_argument("-p","--pixel_size",help="Pixel size for output product",type=float,default=30.0)
  parser.add_argument("-g","--gamma0",action="store_true",help="Make output gamma0 instead of sigma0")
  args = parser.parse_args()

  logFile = "{}_{}_log.txt".format(args.outfile,os.getpid())
  logging.basicConfig(filename=logFile,format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG)
  logging.getLogger().addHandler(logging.StreamHandler())
  logging.info("Starting run")

  geocode_sentinel(args.infile,args.outfile,height=args.terrain_height,pixel_size=args.pixel_size,gamma0_flag=args.gamma0)

