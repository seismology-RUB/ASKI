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
from netCDF4 import Dataset  # https://github.com/Unidata/netcdf4-python
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
            vtk_fid.write("GMT Shorelines, projected / cut to inversion grid\n")
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
            vtk_fid.write("GMT Shorelines, projected / cut to inversion grid\n")
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
    parser = argparse.ArgumentParser(description='create shorelines vtk file with data from gmt file, projected on current inversion grid surface')

    # add options
    parser.add_argument('-recr','--recreate',action='store_true',required=False, 
                        help="if set, then existing inversion grids will be recreated and not just read")
    # add positional arguments
    parser.add_argument('filename_gmt',action='store',
                        help="NetCDF file name containing gmt shore lines data")
    parser.add_argument('outfile_base',action='store',
                        help="absolute output file base (will be concatenated by '.vtk' for vtk output)")
    parser.add_argument('main_parfile',action='store',
                        help="main parameter file of inversion")
    args = parser.parse_args()

    #NetCDF_file = '/home/florian/work/gmt_shoreline_data/binned_GSHHS_c.cdf'
    #NetCDF_file = '/home/florian/work/gmt_shoreline_data/binned_GSHHS_l.cdf'
    #NetCDF_file = '/home/florian/work/gmt_shoreline_data/binned_GSHHS_i.cdf'
    #NetCDF_file = '/home/florian/work/gmt_shoreline_data/binned_GSHHS_h.cdf'
    #NetCDF_file = '/home/florian/work/gmt_shoreline_data/binned_GSHHS_f.cdf'
    NetCDF_file = args.filename_gmt

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

    print "\n+++ opening NetCDF file '"+NetCDF_file+"'"

    NetCDF_fid = Dataset(NetCDF_file, 'r')

    print "\n+++ reading in content of NetCDF file"

    bin_size = NetCDF_fid.variables["Bin_size_in_minutes"][:]
    bin_nx = NetCDF_fid.variables["N_bins_in_360_longitude_range"][:]
    bin_ny = NetCDF_fid.variables["N_bins_in_180_degree_latitude_range"][:]
    n_bin = NetCDF_fid.variables["N_bins_in_file"][:]
    n_seg = NetCDF_fid.variables["N_segments_in_file"][:]
    n_pt = NetCDF_fid.variables["N_points_in_file"][:]
    bin_firstseg = NetCDF_fid.variables["Id_of_first_segment_in_a_bin"][:]
    bin_info = NetCDF_fid.variables["Embedded_node_levels_in_a_bin"][:]
    bin_nseg = NetCDF_fid.variables["N_segments_in_a_bin"][:]
    seg_info = NetCDF_fid.variables["Embedded_npts_levels_exit_entry_for_a_segment"][:]
    GSHHS_area = NetCDF_fid.variables["Ten_times_the_km_squared_area_of_polygons"][:]
    seg_start = NetCDF_fid.variables["Id_of_first_point_in_a_segment"][:]
    pt_dx = NetCDF_fid.variables["Relative_longitude_from_SW_corner_of_bin"][:]
    pt_dy = NetCDF_fid.variables["Relative_latitude_from_SW_corner_of_bin"][:]
    GSHHS_areafrac = NetCDF_fid.variables["Micro_fraction_of_full_resolution_area"][:]
    n_poly = NetCDF_fid.variables["N_polygons_in_file"][:]
    n_node = NetCDF_fid.variables["N_nodes_in_file"][:]
    GSHHS_parent = NetCDF_fid.variables["Id_of_parent_polygons"][:]
    GSHHS_node = NetCDF_fid.variables["Id_of_node_polygons"][:]
    seg_GSHHS_ID = NetCDF_fid.variables["Id_of_GSHHS_ID"][:]
    scale = (bin_size[0]/60.0)/65535.0
    bsize = bin_size[0]/60.0

    print "  bin size in deg = "+str(bsize)
    print "  N bins in x = "+str(bin_nx[0])+" ; N bins in y = "+str(bin_ny[0])+" ; number of bins = "+str(n_bin[0])
    print "  number of segments = "+str(n_seg[0])+" ; number of points = "+str(n_pt[0])

    # select which bins to loop on (i.e. discard bins that do not need to be checked)
    # until implementing an actual selection, select all bins
    idx_bins = range(n_bin[0])

    print "\n+++ looping on "+str(len(idx_bins))+" bins"

    lat_vtk = []
    lon_vtk = []
    lines_vtk = []
    np_vtk = 0 # initialize global counter of vtk points
    for ib in idx_bins:
        #print "  bin index in file: "+str(ib)
        lonsw = numpy.mod(ib,bin_nx[0])*bin_size[0]/60.0 # lower left corner of bin
        latsw = 90.0-((ib/bin_nx[0])+1)*bin_size[0]/60.0

        #print "    looping on "+str(bin_nseg[ib])+" segments contained in this bin"
        for iseg in range(bin_firstseg[ib],bin_firstseg[ib]+bin_nseg[ib]):
            lat_seg = []
            lon_seg = []
            ichunk_now = -1
            #print "seg_info = "+str(seg_info[iseg])+" ; binary representation = "+numpy.binary_repr(seg_info[iseg])
            np_seg = numpy.right_shift(seg_info[iseg],9)
            #print "np_seg = "+str(np_seg)+" ; binary representation = "+numpy.binary_repr(np_seg)
            #print pt_dx(iseg)
            dx = pt_dx[seg_start[iseg]:seg_start[iseg]+np_seg]
            dy = pt_dy[seg_start[iseg]:seg_start[iseg]+np_seg]

            for ip_seg in range(np_seg):
                if dx[ip_seg]<0:
                    lon = lonsw+(65536+dx[ip_seg])*scale               # dx is unsigned short, so shift "negative" values of dx by 2^16
                else:
                    lon = lonsw+dx[ip_seg]*scale
                if dy[ip_seg]<0:
                    lat = latsw+(65536+dy[ip_seg])*scale               # dy is unsigned short, so shift "negative" values of dx by 2^16
                else:
                    lat = latsw+dy[ip_seg]*scale

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
                            # create a new empty segment in which the current point will be stored
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
                        # and create a new empty vtk segment for rest of the gmt segment (possibly gmt segment enters inversion grid again)
                        lat_vtk += lat_seg
                        lon_vtk += lon_seg
                        lines_vtk.append( range(np_vtk,np_vtk+len(lat_seg)) ) # vtk line is defined by list of vtk point indices (starting off by 0)
                        np_vtk += len(lat_seg) # increase global counter of vtk points
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
                ####xyz_vtk += xyz_seg
                lines_vtk.append( range(np_vtk,np_vtk+len(lat_seg)) ) # vtk line is defined by list of vtk point indices (starting off by 0)
                np_vtk += len(lat_seg) # increase global counter of vtk points

    # print "\n+++ create plot axis"
    # plt.axis([0,360,-90,90])
    # print "\n+++ loop on "+str(len(x))+" segments to plot"
    # for xseg,yseg in zip(x,y):
    #     plt.plot(xseg,yseg,'k')
    # print "\n+++ show plot"
    # plt.show()

    print "\n+++ closing NetCDF file"

    NetCDF_fid.close()

    if(np_vtk == 0):
        print "\n+++ NO shore lines were found (no sensible (partial) segments inside the inversion grid domain) -> NO vtk file will be produced"
        exit_now()

    print "\n+++ transforming "+str(np_vtk)+" lat/lon pairs extracted from  NetCDF file to coordinates in vtk projection '"+invgrid_vtk_projection+"'"
    print "  "+str(len(lines_vtk))+" vtk polygons were found"
    print "  (possibly some shore line segments were split up/truncated in order to fit inside the inversion grid domain and not to intersect chunk boundaries in an unwanted manner)"

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
