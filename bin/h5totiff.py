#!/usr/bin/env python3
"""
v1.2.1 20210107 Yu Morishita, GSI

========
Overview
========
This script outputs a float32 file of cumulative displacement from cum*.h5.

=====
Usage
=====
h5totiff.py -d yyyymmdd [-i infile] [-o outfile] [-m yyyymmdd] [-r x1:x2/y1:y2]
     [--ref_geo lon1/lon2/lat1/lat2] [--mask maskfile] [--png] 

 -d  Date to be output
 -i  Path to input cum file (Default: cum_filt.h5)
 -o  Output float32 file (Default: yyyymmdd_yyyymmdd.cum)
 -m  Master (reference) date (Default: first date)
 -r  Reference area (Default: same as info/*ref.txt)
     Note: x1/y1 range 0 to width-1, while x2/y2 range 1 to width
     0 for x2/y2 means all. (i.e., 0:0/0:0 means whole area).
 --ref_geo  Reference area in geographical coordinates.
 --mask  Path to mask file for ref phase calculation (Default: No mask)
 --png   Make png file (Default: Not make png)

"""
#%% Change log
'''
v1.1 20210507 Takayuki Shinoahra, Tokyo Tech
 - Export as Geotiff without world file, export mask data as geotif
v1.0 20210505 Takayuki Shinoahra, Tokyo Tech
 - Original implementation
'''

#%% Import
import getopt
import os
import sys
import re
import time
import numpy as np
import h5py as h5
import SCM
import LiCSBAS_io_lib as io_lib
import LiCSBAS_tools_lib as tools_lib
import LiCSBAS_plot_lib as plot_lib
import imageio
from osgeo import gdal, gdalconst, gdal_array, ogr


class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% Main
def main(argv=None):
   
    #%% Check argv
    if argv == None:
        argv = sys.argv
        
    start = time.time()
    ver="1.0"; date=20210505; author="T. Shinohara"
    print("\n{} ver{} {} {}".format(os.path.basename(argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(argv[0]), ' '.join(argv[1:])), flush=True)


    #%% Set default
    imd_s = []
    cumfile = 'cum_filt.h5'
    outfile = []
    imd_m = []
    refarea = []
    refarea_geo = []
    maskfile = []
    pngflag = False
    cmap = SCM.roma.reversed()


    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hd:i:o:m:r:", ["help", "png", "ref_geo=", "mask="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-d':
                imd_s = a
            elif o == '-i':
                cumfile = a
            elif o == '-o':
                outfile = a
            elif o == '-m':
                imd_m = a
            elif o == '-r':
                refarea = a
            elif o == '--ref_geo':
                refarea_geo = a
            elif o == '--mask':
                maskfile = a
            elif o == '--png':
                pngflag = True

        # if not imd_s:
        #     raise Usage('No date given, -d is not optional!')
        if not os.path.exists(cumfile):
            raise Usage('No {} exists! Use -i option.'.format(cumfile))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2


    #%% Read info
    ### Read cumfile
    cumh5 = h5.File(cumfile,'r')
    imdates = cumh5['imdates'][()].astype(str).tolist()
    cum = cumh5['cum']
    n_im, length, width = cum.shape



    ### Reference area
    if refarea:
        if not tools_lib.read_range(refarea, width, length):
            print('\nERROR in {}\n'.format(refarea), file=sys.stderr)
            return 2
        else:
            refx1, refx2, refy1, refy2 = tools_lib.read_range(refarea, width, length)
    elif refarea_geo:
        lat1 = float(cumh5['corner_lat'][()])
        lon1 = float(cumh5['corner_lon'][()])
        dlat = float(cumh5['post_lat'][()])
        dlon = float(cumh5['post_lon'][()])
        if not tools_lib.read_range_geo(refarea_geo, width, length, lat1, dlat, lon1, dlon):
            print('\nERROR in {}\n'.format(refarea_geo), file=sys.stderr)
            return 2
        else:
            refx1, refx2, refy1, refy2 = tools_lib.read_range_geo(refarea_geo, width, length, lat1, dlat, lon1, dlon)
    else:
        refarea = cumh5['refarea'][()]
        refx1, refx2, refy1, refy2 = [int(s) for s in re.split('[:/]', refarea)]
    
    ### Master (reference) date
    if not imd_m:
        imd_m = imdates[0]
        
    ### mask
    if maskfile:
        mask = io_lib.read_img(maskfile, length, width)
        # mask[mask==0] = np.nan
    else:
        mask = np.ones((length, width), dtype=np.float32)
        ix_m = imdates.index(imd_m)
        mask[np.isnan(cum[ix_m, :, :])] = np.nan
    
    ### save mask
    # # save as geotiff
    dlat = -0.0009999992325901985
    dlon = 0.0009999992325901985
    lat_n_g = float(cumh5['corner_lat'][()]) #grid reg
    lon_w_g = float(cumh5['corner_lon'][()]) #grid reg∂

    ## Grid registration to pixel registration by shifing half pixel
    lat_n_p = lat_n_g - dlat/2
    lon_w_p = lon_w_g - dlon/2
    compress_option = ['COMPRESS=DEFLATE', 'PREDICTOR=3']
    nodata = np.nan
    io_lib.make_geotiff(mask, lat_n_p, lon_w_p, dlat, dlon, "./mask.tif", compress_option, nodata)

    ### 全部のスレイブ画像を処理する
    for i in range(len(imdates)):
        imd_s=imdates[i]
        # print(img_s)
        ix_s = imdates.index(imd_s)
        ix_m = imdates.index(imd_m)

        ### Check date
        if not imd_s in imdates:
            print('\nERROR: No date of {} exist in {}!'.format(imd_s, cumfile), file=sys.stderr)
            return 2
        if not imd_m in imdates:
            print('\nERROR: No date of {} exist in {}!'.format(imd_m, cumfile), file=sys.stderr)
            return 2


        outfile = '{}_{}.cum'.format(imd_m, imd_s)
        print(outfile)
        
        #%% Make flt
        cum_s = cum[ix_s, :, :]
        cum_m = cum[ix_m, :, :]

        cum_dif = cum_s-cum_m
        cum_dif = cum_dif-np.nanmean(cum_dif[refy1:refy2, refx1:refx2])
        cum_dif = cum_dif*mask
            
        cum_dif.tofile(outfile)

        io_lib.make_geotiff(cum_dif, lat_n_p, lon_w_p, dlat, dlon, outfile.replace("cum", "tif") , compress_option, nodata)

        #%% Make png if specified
        if pngflag:
            pngfile = outfile+'.png'
            title = '{} (Ref X/Y {}:{}/{}:{})'.format(outfile, refx1, refx2, refy1, refy2)
            plot_lib.make_im_png(cum_dif, pngfile, cmap, title)




    #%% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}\n'.format(outfile), flush=True)


#%% main
if __name__ == "__main__":
    sys.exit(main())
