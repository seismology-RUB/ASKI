#!/usr/bin/env python
#
#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
import argparse
from inputParameter import inputParameter
from os import path as os_path
import numpy
from create_shore_lines_f2py import shore_f2py
from sys import exit as sys_exit
#
#--------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#--------------------------------------------------------------------------------
def exit_now():
    shore_f2py.deallocate_invgrid()
    print "\n+++ Good Bye\n"
    sys_exit()
#
#--------------------------------------------------------------------------------
#
# THE FOLLOWING DOCUMENTATION IS TAKEN FROM 
# https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/readme.txt
#
#
# struct GSHHG {  /* Global Self-consistent Hierarchical High-resolution Shorelines */
#         int id;         /* Unique polygon id number, starting at 0 */
#         int n;          /* Number of points in this polygon */
#         int flag;       /* = level + version << 8 + greenwich << 16 + source << 24 + river << 25 */
#         /* flag contains 5 items, as follows:
#          * low byte:    level = flag & 255: Values: 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
#          * 2nd byte:    version = (flag >> 8) & 255: Values: Should be 12 for GSHHG release 12 (i.e., version 2.2)
#          * 3rd byte:    greenwich = (flag >> 16) & 1: Values: Greenwich is 1 if Greenwich is crossed
#          * 4th byte:    source = (flag >> 24) & 1: Values: 0 = CIA WDBII, 1 = WVS
#          * 4th byte:    river = (flag >> 25) & 1: Values: 0 = not set, 1 = river-lake and level = 2
#          */
#         int west, east, south, north;   /* min/max extent in micro-degrees */
#         int area;       /* Area of polygon in 1/10 km^2 */
#         int area_full;  /* Area of original full-resolution polygon in 1/10 km^2 */
#         int container;  /* Id of container polygon that encloses this polygon (-1 if none) */
#         int ancestor;   /* Id of ancestor polygon in the full resolution set that was the source of this polygon (-1 if none) */
# };
#
# Following each header structure is n structures of coordinates:
#
# struct GSHHG_POINT {	/* Each lon, lat pair is stored in micro-degrees in 4-byte signed integer format */
# 	int32_t x;
# 	int32_t y;
# };
#
# Some useful information:
#
# A) To avoid headaches the binary files were written to be big-endian.
#    If you use the GMT supplement gshhg it will check for endian-ness and if needed will
#    byte swab the data automatically. If not then you will need to deal with this yourself.
#
# B) In addition to GSHHS we also distribute the files with political boundaries and
#    river lines.  These derive from the WDBII data set.
#
# C) As to the best of our knowledge, the GSHHG data are geodetic longitude, latitude
#    locations on the WGS-84 ellipsoid.  This is certainly true of the WVS data (the coastlines).
#    Lakes, riverlakes (and river lines and political borders) came from the WDBII data set
#    which may have been on WGS072.  The difference in ellipsoid is way less then the data
#    uncertainties.  Offsets have been noted between GSHHG and modern GPS positions.
#
# D) Originally, the gshhs_dp tool was used on the full resolution data to produce the lower
#    resolution versions.  However, the Douglas-Peucker algorithm often produce polygons with
#    self-intersections as well as create segments that intersect other polygons.  These problems
#    have been corrected in the GSHHG lower resolutions over the years.  If you use gshhs_dp to
#    generate your own lower-resolution data set you should expect these problems.
#
# E) The shapefiles release was made by formatting the GSHHG data using the extended GMT/GIS
#    metadata understood by OGR, then using ogr2ogr to build the shapefiles.  Each resolution
#    is stored in its own subdirectory (e.g., f, h, i, l, c) and each level (1-4) appears in
#    its own shapefile.  Thus, GSHHS_h_L3.shp contains islands in lakes for the high res
#    data. Because of GIS limitations some polygons that straddle the Dateline (including
#    Antarctica) have been split into two parts (east and west).
#
# F) The netcdf-formatted coastlines distributed with GMT derives directly from GSHHG; however
#    the polygons have been broken into segments within tiles.  These files are not meant
#    to be used by users other than via GMT tools (pscoast, grdlandmask, etc).
def read_shore_line_file(filename):
    polygons = []
    npol = 0

    with open(filename,'rb') as f:
        # read polygon header, consisting of 11 4-byte integers
        header_byte_string = f.read(44)

        while header_byte_string:
            # check if header is shorter than 44 bytes
            if len(header_byte_string) != 44:
                print ("### ERROR : header of "+str(npol)+"'th polygon has length of "+str(len(header_byte_string))+
                       " bytes; expecting 44 bytes, i.e. file ended in the middle of a polygon header => unexpected file format!\n")
                return []

            # header consists of 11 4-byte signed integers (big-endian)
            header_int = numpy.fromstring(header_byte_string,dtype=numpy.dtype(">i4"))

            if npol != header_int[0]:
                # first entry of header does not coincide with polygon ID = index starting from 0
                print ("### ERROR : the header of the "+str(npol)+"'th polygon says: polygon ID = "+str(header_int[0])+
                       "; expecting value "+str(npol)+" => unexpected file format!\n")
                return []

            # second entry of header is number of points of polygon
            npoint = header_int[1]
            if npoint <=0:
                print ("### ERROR : the header of the "+str(npol)+"'th polygon says: number of points in polygon = "+str(npoint)+
                       "; expecting strictly positive value => unexpected file format!\n")
                return []

            # read byte string of polygon lon/lat pairs
            polygon_byte_string = f.read(4*2*npoint)
            if len(polygon_byte_string) != 4*2*npoint:
                print ("### ERROR : reading "+str(npol)+"'th polygon's points: could only read "+str(len(polygon_byte_string))+
                       " bytes from file, but expecting "+str(npoint)+" points (as header says), i.e. 8*"+str(npoint)+" = "+str(8*npoint)+
                       " bytes; file ended in the middle of this polygon => unexpected file format!\n")
                return []
            # polygon lon/lat pairs consist of signed integers (big-endian) and represent micro degrees
            polygon_int = numpy.fromstring(polygon_byte_string,dtype=numpy.dtype(">i4"))
            lon_lat_list = [(1.e-6*polygon_int[2*i],1.e-6*polygon_int[2*i+1]) for i in range(npoint)]
            polygons.append(lon_lat_list)

            npol += 1

            # try to read next polygon header, consisting of 11 4-byte integers
            # if the file ends here, header_byte_string will be '' and loop will terminate
            header_byte_string = f.read(44)

    return polygons
#
#--------------------------------------------------------------------------------
#
def split_vtk_segment(ichunk_previous,ichunk_now,invgrid_type,invgrid_vtk_projection):
    #print "ichunk_previous,ichunk_now,invgrid_type,invgrid_vtk_projection = ",ichunk_previous,ichunk_now,invgrid_type,invgrid_vtk_projection
    if ichunk_previous <= 0 or ichunk_now <= 0:
        #print "## at least one of ichunk_previous,ichunk_now are -1"
        return False
    if ichunk_previous == ichunk_now:
        #print "## ichunk_previous and ichunk_now are equal"
        return False

    # in the following, assume that ichunk_previous != ichunk_now and that both are > 0

    if invgrid_type == 'chunksInversionGrid':
        if not "FLAT" in invgrid_vtk_projection:
            return False

        #
        # DISTRIBUTION OF CHUNKS IN chunksInversionGrid (LOCAL_FLAT projection here):
        #
        #     +---+
        #     | 3 |
        # +---+---+---+---+
        # | 2 | 1 | 5 | 4 |
        # +---+---+---+---+
        #     | 6 |
        #     +---+
        #
        if ichunk_previous == 1:
            split_seg = ichunk_now not in [2,3,5,6]
            #print "## chunksInversionGrid: ichunk_previous = 1, ichunk_now = "+str(ichunk_now)+", split_seg = ",split_seg
        elif ichunk_previous in [2,3,6]:
            split_seg = ichunk_now not in [1]
            #print "## chunksInversionGrid: ichunk_previous = "+str(ichunk_previous)+", ichunk_now = "+str(ichunk_now)+", split_seg = ",split_seg
        elif ichunk_previous == 4:
            split_seg = ichunk_now not in [5]
        elif ichunk_previous == 5:
            split_seg = ichunk_now not in [1,4]
        return split_seg
    else:
        return False
#
#--------------------------------------------------------------------------------
#
def vtk_file_exists(filebase):
    valid_name = filebase+'.vtk'
    icount = 0
    while os_path.exists(valid_name):
        icount += 1
        valid_name = filebase+'-%i.vtk'%icount
    file_existed = icount > 0
    return file_existed,valid_name
#
#--------------------------------------------------------------------------------
#
def write_vtk_file(filebase,vtk_format,points,lines):
    file_exists,filename = vtk_file_exists(filebase)
    if file_exists:
        print "  WARNING: the output file '"+filebase+".vtk' already exists; adding a numeral extension"
        print "  writing "+vtk_format+" vtk file '"+filename+"'"
        pass
    else:
        print "  writing "+vtk_format+" vtk file '"+filename+"'"

    # ASCII vtk file
    if vtk_format == 'ASCII':
        with open(filename,'w') as vtk_fid:
            # HEADER
            vtk_fid.write("# vtk DataFile Version 3.1\n")
            vtk_fid.write("GSHHS Shorelines, projected / cut to inversion grid\n")
            vtk_fid.write("ASCII\n")
            vtk_fid.write("DATASET POLYDATA\n")

            # POINTS
            vtk_fid.write("POINTS %i float\n"%len(points))
            for point in points:
                vtk_fid.write("%f %f %f\n"%point)

            # LINES
            number_of_lines_values = sum([1+len(line) for line in lines])
            vtk_fid.write("LINES %i %i\n"%(len(lines),number_of_lines_values))
            for line in lines:
                vtk_fid.write(str(len(line))+" "+" ".join([str(point_index) for point_index in line])+"\n")
    else:
        with open(filename,'wb') as vtk_fid:
            # HEADER
            vtk_fid.write("# vtk DataFile Version 3.1\n")
            vtk_fid.write("GSHHS Shorelines, projected / cut to inversion grid\n")
            vtk_fid.write("BINARY\n")
            vtk_fid.write("DATASET POLYDATA\n")

            # POINTS
            vtk_fid.write("POINTS %i float\n"%len(points))
            for point in points:
                # create numpy array of BigEndian 4-byte floats (binary VTK legacy format must be written in BigEndian)a
                np_point = numpy.array(point,dtype=numpy.dtype(">f4"))
                vtk_fid.write(np_point.tobytes())
            vtk_fid.write("\n")

            # LINES
            number_of_lines_values = sum([1+len(line) for line in lines])
            vtk_fid.write("LINES %i %i\n"%(len(lines),number_of_lines_values))
            for line in lines:
                # create numpy array of BigEndian 4-byte integers (binary VTK legacy format must be written in BigEndian)a
                np_line = numpy.array([len(line)]+line,dtype=numpy.dtype(">i4"))
                vtk_fid.write(np_line.tobytes())
            vtk_fid.write("\n")
#--------------------------------------------------------------------------------
# MAIN PROGRAM
#--------------------------------------------------------------------------------
def main():
    # initialize argument parser 
    parser = argparse.ArgumentParser(description='create shorelines vtk file with data from GSHHS file, projected on current inversion grid surface')

    # add options
    parser.add_argument('-recr','--recreate',action='store_true',required=False, 
                        help="if set, then existing inversion grids will be recreated and not just read")
    # add positional arguments
    parser.add_argument('GSHHS_filename',action='store',
                        help="native binary GSHHS shore line file")
    parser.add_argument('outfile_base',action='store',
                        help="absolute output file base (will be concatenated by '.vtk' for vtk output)")
    parser.add_argument('main_parfile',action='store',
                        help="main parameter file of inversion")
    args = parser.parse_args()

    native_binary_file = args.GSHHS_filename

    vtk_outfile_base = args.outfile_base

    main_parfile = args.main_parfile

    invgrid_recreate = args.recreate

    print "\n+++ welcome to shore line creation for ASKI!"

    print "\n+++ retrieving all information from main parameter file '"+main_parfile+"'"

    mparam = inputParameter(main_parfile)
    noKeys = mparam.keysNotPresent(['MAIN_PATH_INVERSION','ITERATION_STEP_PATH','CURRENT_ITERATION_STEP',
                                    'PARFILE_ITERATION_STEP','DEFAULT_VTK_FILE_FORMAT'])
    if len(noKeys) > 0:
        print "### ERROR : the following keywords are required in main parameter file '"+main_parfile+"':"
        print "### "+',  '.join(noKeys)
        exit_now()

    vtk_format = mparam.sval('DEFAULT_VTK_FILE_FORMAT')
    if not vtk_format in ['ASCII','BINARY']:
        print "### ERROR : DEFAULT_VTK_FILE_FORMAT = "+vtk_format+" in main parameter file. Must be either ASCII or BINARY"
        exit_now()

    iter_path = os_path.join(mparam.sval('MAIN_PATH_INVERSION'),mparam.sval('ITERATION_STEP_PATH')+
                                '%3.3i/'%mparam.ival('CURRENT_ITERATION_STEP'))
    iter_parfile = os_path.join(iter_path,mparam.sval('PARFILE_ITERATION_STEP'))

    print "  iteration step parfile = '"+iter_parfile+"'"
    print "  default vtk file format (will be used here) = "+vtk_format

    print "\n+++ retrieving all information from iteration step parameter file '"+iter_parfile+"'"

    iparam = inputParameter(iter_parfile)
    noKeys = iparam.keysNotPresent(['TYPE_INVERSION_GRID','PARFILE_INVERSION_GRID'])
    if len(noKeys) > 0:
        print "### ERROR : the following keywords are required in iteration-step  parameter file '"+iter_parfile+"':"
        print "### "+',  '.join(noKeys)
        exit_now()

    invgrid_type = iparam.sval('TYPE_INVERSION_GRID')
    invgrid_parfile = os_path.join(iter_path,iparam.sval('PARFILE_INVERSION_GRID'))
    invgrid_path = iter_path

    print "  inversion grid type = "+invgrid_type
    print "  inversion grid parameter file = "+invgrid_parfile

    print "\n+++ retrieving all information from iteration step parameter file '"+iter_parfile+"'"

    igparam = inputParameter(invgrid_parfile)
    noKeys = igparam.keysNotPresent(['VTK_PROJECTION'])
    if len(noKeys) > 0:
        print "### ERROR : the following keywords are required in inversion grid  parameter file '"+invgrid_parfile+"':"
        print "### "+',  '.join(noKeys)
        exit_now()

    invgrid_vtk_projection = igparam.sval('VTK_PROJECTION')

    print "  vtk projection of inversion grid = "+invgrid_vtk_projection

    if invgrid_recreate:
        print "\n+++ setting up the inversion grid, recreating it if existend (as indicated by option -recr)"
    else:
        print "\n+++ setting up the inversion grid, will only read it in, if existend (option -recr is NOT set)"
    status = shore_f2py.set_up_invgrid(invgrid_type,invgrid_parfile,invgrid_path,invgrid_recreate)
    if status != 0:
        print "### ERROR : inversion grid could not be set up"
        exit_now()

    print "\n+++ reading native-binary GSHHS file '"+native_binary_file+"'"

    try:
        polygons = read_shore_line_file(native_binary_file)
    except:
        print "### ERROR : could not open the file\n"
        shore_f2py.deallocate_invgrid()
        raise
    print "  found "+str(len(polygons))+" polygons"

    #print "  bin size in deg = "+str(bsize)
    #print "  N bins in x = "+str(bin_nx[0])+" ; N bins in y = "+str(bin_ny[0])+" ; number of bins = "+str(n_bin[0])
    #print "  number of segments = "+str(n_seg[0])+" ; number of points = "+str(n_pt[0])


    lat_vtk = []
    lon_vtk = []
    ####xyz_vtk = []
    lines_vtk = []
    np_vtk = 0 # initialize global counter of vtk points
    for polygon in polygons:
        lat_seg = []
        lon_seg = []
        ichunk_now = -1
        for lon,lat in polygon:

            # test if point is inside inversion grid volume
            ichunk_previous = ichunk_now
            status,ichunk_now = shore_f2py.point_is_inside_inversion_grid(lat,lon)

            # YES, point is inside the inverison grid volume:
            if status == 0:
                # if this is not the first point of the vtk segment inside the inversion grid, check if it
                # crosses a chunk border which is discontinuous in the current inversion grid's vtk projection
                # so that the vtk segment needs to be split here
                if(len(lat_seg)>=1):
                    # compare ichunk_now with ichunk_previous (of last point in lat/lon_seg) and decide whether the segment has to be split here (dependent on ichunk values and invgrid type and vtk projection)
                    if split_vtk_segment(ichunk_previous,ichunk_now,invgrid_type,invgrid_vtk_projection):
                        # if vtk segment needs to be splitted here, check whether there was a sensible vtk line found so far (at least 2 points)
                        if(len(lat_seg)>=2):
                            # if there was a sensible vtk line found so far, memorize it
                            lat_vtk += lat_seg
                            lon_vtk += lon_seg
                            lines_vtk.append( range(np_vtk,np_vtk+len(lat_seg)) ) # vtk line is defined by list of vtk point indices (starting off by 0)
                            np_vtk += len(lat_seg) # increase global counter of vtk points
                        else:
                            # if there is only one point found so far, throw it away, i.e. do not memorize anything
                            pass
                        # create a new empty segment in which the current point will be stored (as first point)
                        lat_seg = []
                        lon_seg = []

                # store point in the current vtk segment (existing one or the one just newly created)
                lat_seg.append(lat)
                lon_seg.append(lon)

            # NO, point is not inside the inversion grid volume:
            else:
                # check whether there was a sensible vtk line found so far (at least 2 points)
                if(len(lat_seg)>=2):
                    # if there was a sensible vtk line found so far, memorize it
                    lat_vtk += lat_seg
                    lon_vtk += lon_seg
                    lines_vtk.append( range(np_vtk,np_vtk+len(lat_seg)) ) # vtk line is defined by list of vtk point indices (starting off by 0)
                    np_vtk += len(lat_seg) # increase global counter of vtk points

                # create a new empty vtk segment for rest of the GSHHS segment (possibly GSHHS segment enters inversion grid again)
                lat_seg = []
                lon_seg = []

                # in any case, make the current chunk index invalid (this point is not inside the inversion grid)
                ichunk_now = -1

        # after looping on all points of the segment, memorize (the tail part of) it for vtk output
        # ONLY IN CASE there are at least 2 points in it.
        # this could be only the tail part of the segment, if some segment points were outside the inversion grid
        if(len(lat_seg)>=2):
            lat_vtk += lat_seg
            lon_vtk += lon_seg
            lines_vtk.append( range(np_vtk,np_vtk+len(lat_seg)) ) # vtk line is defined by list of vtk point indices (starting off by 0)
            np_vtk += len(lat_seg) # increase global counter of vtk points

    # print "\n+++ create plot axis"
    # plt.axis([0,360,-90,90])
    # print "\n+++ loop on "+str(len(x))+" segments to plot"
    # for xseg,yseg in zip(x,y):
    #     plt.plot(xseg,yseg,'k')
    # print "\n+++ show plot"
    # plt.show()


    if(np_vtk == 0):
        print "\n+++ NO shore lines were found (no sensible (partial) segments inside the inversion grid domain) -> NO vtk file will be produced"
        exit_now()

    print "\n+++ transforming "+str(np_vtk)+" lat/lon pairs extracted from native binary file to coordinates in vtk projection '"+invgrid_vtk_projection+"'"
    print "  "+str(len(lines_vtk))+" vtk line segments were generated"
    print "  (possibly some original shore line polygons were split up/truncated/removed in order to fit inside the inversion grid domain and not to intersect chunk boundaries in an unwanted manner)"

    status,x_vtk,y_vtk,z_vtk = shore_f2py.transform_to_vtk_projection(lat_vtk,lon_vtk)
    if status != 0:
        print "### ERROR : lat/lon coordinates could not be transformed to vtk projection"
        exit_now()
    xyz_vtk = zip(x_vtk,y_vtk,z_vtk)
    #xyz_vtk = zip(lon_vtk,lat_vtk,[0.0 for i in range(np_vtk)])

    print "\n+++ writing vtk output file"

    write_vtk_file(vtk_outfile_base,vtk_format,xyz_vtk,lines_vtk)

    exit_now()


if __name__=="__main__":
    main()
