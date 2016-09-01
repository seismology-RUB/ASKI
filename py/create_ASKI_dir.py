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
from sys import argv as sys_argv
from sys import exit as sys_exit
from os import path as os_path
from os import mkdir as os_mkdir
from os import makedirs as os_makedirs
from inputParameter import inputParameter
#
#############################################################
# define names of files and paths here to be used to create
# iter_parfile
# if you want to use a different naming, change here
FILE_WAVEFIELD_POINTS = ''
FILE_KERNEL_REFERENCE_MODEL = ''

PATH_OUTPUT_FILES = 'output_files/'

PATH_KERNEL_DISPLACEMENTS = 'kernel_displacements/'
PATH_KERNEL_GREEN_TENSORS = 'kernel_green_tensors/'
PATH_SENSITIVITY_KERNELS = 'sensitivity_kernels/'
PATH_SYNTHETIC_DATA = 'synthetic_data/'
PATH_KERNEL_REFERENCE_MODELS = 'kernel_reference_models/'

FILE_INTEGRATION_WEIGHTS = 'integration_weights'
FILEBASE_BASIC_STATS = PATH_OUTPUT_FILES+'stats'
#
#############################################################
def exit_now():
    print "###"
    print "Good Bye"
    print "###"
    print ""
    sys_exit()
#############################################################
#
def create_iter_dir(iter_path,indx):
    iter_parfile_str = """#------------------------------------------------------------------------------
#  PARAMETER FILE FOR INVERSION "{0}", ITERATION STEP {1}
#------------------------------------------------------------------------------
#
#  HOWTO USE THE PARAMETER FILE:
#
# comment lines, STARTING with "#" are ignored
# empty lines and lines containing blanks only are ignored
# lines not containing "=" are ignored
#
# specify one parameter per line, valid lines having the form "keyword = value"
# (blanks leading or following the keyword, the "=", or the value are ignored)
#
# in a non comment line, all characters in front of "=" (without leading and appending blanks)
# are interpreted as the keyword, allowing for blank characters within the keyword (e.g. for lines  "   key word = value   ",
# the string "key word" is used as the keyword)
# all characters behind "=" (without leading and appending blanks) are interpreted as the value string from 
# which the value is read, which in particular means that "#" comments at the end of a line (such as 
# " keyword = value  # comment ") are NOT allowed (in the latter case, the string "value  # comment" would be used
# to read a value from)
#
# by ASKI convention, specify PATHS (i.e. directory names, which will be 
# concatenated with a filename of a file in that directory) always ending on '/'
# specify FILENAMES always WITHOUT leading '/'
#------------------------------------------------------------------------------

#---------------------------------------------------------------------------
#  PATH SPECIFIC MODE
#---------------------------------------------------------------------------

#  set the following flag to .true. if you intend to use the ASKI functionality of 
#  path specific kernel reference models (usually only for the very first iteration
#  and a 1D method)
#  if you want to use one global kernel reference model (for all paths, i.e. kernels)
#  set this flag to .false.
#
  USE_PATH_SPECIFIC_MODELS = .false.

#---------------------------------------------------------------------------
#  FREQUENCY SUBSET
#---------------------------------------------------------------------------

#  frequency discretization of this iteration step, must be subset of global frequency discretization 
#  for this inversion defined in main parameter file
#  define subset by ITERATION_STEP_NUMBER_OF_FREQ <= MEASURED_DATA_NUMBER_OF_FREQ
#  and a vector ITERATION_STEP_INDEX_OF_FREQ of frequency indices (of length ITERATION_STEP_NUMBER_OF_FREQ)
#  for which all indices contained in vector ITERATION_STEP_INDEX_OF_FREQ must also be contained in global
#  vector MEASURED_DATA_INDEX_OF_FREQ
#  all indices here are assumed in accordance with the global frequency step MEASURED_DATA_FREQUENCY_STEP
#
  ITERATION_STEP_NUMBER_OF_FREQ = 
  ITERATION_STEP_INDEX_OF_FREQ = 


#---------------------------------------------------------------------------
#  INVERSION GRID
#---------------------------------------------------------------------------

#  Type of inversion grid:
#    chunksInversionGrid
#      ASKI internal, method independent spherical inverison grid using (one or several) chunks in 
#      the concept of a cubed sphere and variable cell refinement
#
#    schunkInversionGrid
#      ASKI internal and method independent spherical simple chunk inverison grid using one chunks in 
#      the concept of a cubed sphere which is assumed rather small in lateral width (as certain
#      assumptions are made on the inversion grid cells not to be distorted too much)
#
#    scartInversionGrid
#       ASKI internal and method independent simple Cartesian inversion grid containing of layeres of 
#       regular hexahedra, layer-dependent refinement possible
#
#    ecartInversionGrid
#       method independent, external Cartesian inversion grid provided e.g. by CUBIT which may contain tetrahedra, 
#       hexahedra are not fully supported yet
#
#    specfem3dInversionGrid
#       method dependent inversion grid for METHOD = SPECFEM3D 
#       uses SPECFEM elements as inversion grid:
#       use spectral elements as the inversion grid cells and read in the geometry including information     
#       on neighbour cells, using ALL GLL points in an element as wavefield points, reading in their jacobian
#
  TYPE_INVERSION_GRID = 


#  specific parameter file used to create an inversion grid of type TYPE_INVERSION_GRID, 
#  relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  PARFILE_INVERSION_GRID = 

#---------------------------------------------------------------------------
#  INTEGRATION WEIGHTS
#---------------------------------------------------------------------------

#  Type of integration weights: 
#    0 = all weights the same, weight = 1/number_of_points_in_box, i.e. no integration, just build average
#    1 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 1
#    2 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 2, i.e. approximation order 3 (?)
#    3 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 3, i.e. approximation order 4 (?)
#    4 = for each cell compute the highest possible order of SDI integration (trying degrees 3,2,1 , in that order). if not even SDI degree 1 is successful, choose weights of type 5
#    5 = average of function values, multiplied with volume of box (i.e. linear integration)
#    6 = external integration weights, to be used along with a suitable inversion grid (e.g. specfem3dInversionGrid)
# 
  TYPE_INTEGRATION_WEIGHTS = 


#  File containing integration weights, relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  FILE_INTEGRATION_WEIGHTS = {9}


#---------------------------------------------------------------------------
#  FILES AND PATHS
#---------------------------------------------------------------------------

#  File containing forward grid points, relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  FILE_WAVEFIELD_POINTS = {2}


#  File containing the reference model for module kernelReferenceModel, relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  FILE_KERNEL_REFERENCE_MODEL = {3}


#  Base filename of vtk stats output files (related to inversion grid, wavefield points, integration weights, 
#  events, stations), relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  FILEBASE_BASIC_STATS = {10}


#  Folder relative to which some sensitivity analysis and inversion programs write their output (relatively small output 
#  like models, coefficients etc., NO wavefields/kernels etc.!), 
#  relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/
#
  PATH_OUTPUT_FILES = {4}


#  Paths, relative to MAIN_PATH_INVERSION/ITERATION_STEP_PATH/ which contain
#  the synthetic wavefields/data and sensitivity kernels
#
#  NAMING CONVENTION OF FILES IN THESE DIRECTORIES, AS EXPECTED BY PROGRAMS:
#  FILE_KERNEL_DISPLACEMENT: will by convention be: kernel_displ_EVENTID (if USE_PATH_SPECIFIC_MODELS = .false.)
#                                                   kernel_displ_EVENTID_STATIONNAME (if USE_PATH_SPECIFIC_MODELS = .true.)
#  FILE_KERNEL_GREEN_TENSOR: will by convention be: kernel_gt_STATIONNAME (if USE_PATH_SPECIFIC_MODELS = .false.)
#                                                   kernel_gt_EVENTID_STATIONNAME (if USE_PATH_SPECIFIC_MODELS = .true.)
#    FILE_KERNEL_GREEN_TENSOR and FILE_KERNEL_DISPLACEMENT will be basenames only.
#    Every method then deals with the naming by its own (could be multiple files or having some extension etc.)
#  FILE_SPECTRAL_KERNEL: will in case of pre-integrated kernels on inversion grid by convention be: spectral_kernel_PARAMETRIZATION_EVENTID_STATIONNAME
#  FILE_SPECTRAL_KERNEL: will in case of plain kernel values on wavefield points by convention be: spectral_kernel_ON-WP_PARAMETRIZATION_EVENTID_STATIONNAME
#  FILE_SYNTHETIC_DATA: must by convention be: synthetics_EVENTID_STATIONNAME_COMP
#  FILE_KERNEL_REFERENCE_MODEL: Only used in case USE_PATH_SPECIFIC_MODELS = .true. ; must by convention be: krm_EVENTID_STATIONNAME
#
  PATH_KERNEL_DISPLACEMENTS = {5}
  PATH_KERNEL_GREEN_TENSORS = {6}
  PATH_SENSITIVITY_KERNELS = {7}
  PATH_SYNTHETIC_DATA = {8}
  PATH_KERNEL_REFERENCE_MODELS = {11}

""".format(inversion_name,indx,FILE_WAVEFIELD_POINTS,FILE_KERNEL_REFERENCE_MODEL,PATH_OUTPUT_FILES,
           PATH_KERNEL_DISPLACEMENTS,PATH_KERNEL_GREEN_TENSORS,PATH_SENSITIVITY_KERNELS,
           PATH_SYNTHETIC_DATA,FILE_INTEGRATION_WEIGHTS,FILEBASE_BASIC_STATS,PATH_KERNEL_REFERENCE_MODELS)
    # create iteration step directory
    print "mkdir '"+iter_path+"'"
    os_mkdir(iter_path)
    # create iter_parfile
    print "create template parameter file '"+os_path.join(iter_path,param.sval('PARFILE_ITERATION_STEP'))+"'"
    open(os_path.join(iter_path,param.sval('PARFILE_ITERATION_STEP')),'w').write(iter_parfile_str)
    # create subdirectories
    print "mkdir '"+os_path.join(iter_path,PATH_OUTPUT_FILES)+"'"
    os_mkdir(os_path.join(iter_path,PATH_OUTPUT_FILES))
    print "mkdir '"+os_path.join(iter_path,PATH_KERNEL_DISPLACEMENTS)+"'"
    os_mkdir(os_path.join(iter_path,PATH_KERNEL_DISPLACEMENTS))
    print "mkdir '"+os_path.join(iter_path,PATH_KERNEL_GREEN_TENSORS)+"'"
    os_mkdir(os_path.join(iter_path,PATH_KERNEL_GREEN_TENSORS))
    print "mkdir '"+os_path.join(iter_path,PATH_SENSITIVITY_KERNELS)+"'"
    os_mkdir(os_path.join(iter_path,PATH_SENSITIVITY_KERNELS))
    print "mkdir '"+os_path.join(iter_path,PATH_SYNTHETIC_DATA)+"'"
    os_mkdir(os_path.join(iter_path,PATH_SYNTHETIC_DATA))
    print "mkdir '"+os_path.join(iter_path,PATH_KERNEL_REFERENCE_MODELS)+"'"
    os_mkdir(os_path.join(iter_path,PATH_KERNEL_REFERENCE_MODELS))
#############################################################
help_str = """
   create_ASKI_dir.py

USAGE: please give 2 arguments:
[1] main parmeter file of the inversion
[2] number of iteration steps

EXAMPLE:
create_ASKI_dir.py ./main_parfile_project1 10
"""
#
#
#############################################################
# START MAIN SCRIPT HERE
#############################################################
#
#
arguments = sys_argv[1:]
#
# check if all required arguments are given
if len(arguments)!=2:
    print help_str
    sys_exit()
parfile = arguments[0]
if not os_path.exists(parfile):
    print "### ERROR: file '"+parfile+"' does not exist"
    print help_str
    sys_exit()
try:
    n_iter_steps = int(arguments[1])
except:
    print "### ERROR: second argument '"+arguments[1]+"' could not be converted to integer"
    print help_str
    sys_exit()
#
print ""
print "###"
print "running script create_ASKI_dir.py now with argument values:"
print "  parameter file = '"+parfile+"'"
print "  number of iteration steps = '"+str(n_iter_steps)+"'"
contin = raw_input('is this correct? continue? [y/n] ')
print "###"
print ""
if not contin.lower().startswith('y'):
    exit_now()
#
param = inputParameter(parfile)
noKeys = param.keysNotPresent(['MAIN_PATH_INVERSION','ITERATION_STEP_PATH','PARFILE_ITERATION_STEP','PARAMETER_CORRELATION_FILE'])
if len(noKeys) > 0:
    print "### ERROR : the following keywords are required in parameter file '"+parfile+"':"
    print "### "+',  '.join(noKeys)
    print ""
    exit_now()
#
main_path = param.sval('MAIN_PATH_INVERSION')
if main_path[-1]=='/':
    main_path = main_path[:-1]
else:
    print ("## WARNING, value '"+main_path+"' of key 'MAIN_PATH_INVERSION' in parfile '"+parfile+
           "' does not end on '/'. Please adjust later (significant for some programs)\n")
inversion_name = os_path.basename(main_path)
#
# now go create directory structure, always checking if paths already exist
print "###"
print ("creating ASKI directory structure with '"+str(n_iter_steps)+
       "' iteration steps in main directory '"+main_path+"'")
print "###"
print ""
#
print "###"
#
# create main path, if it does not exist
if os_path.exists(main_path):
    print "# NOTE, main directory '"+main_path+"' already exists."
else:
    print "mkdir '"+main_path+"'"
    os_makedirs(main_path)
#
# create empty model parameter correlation file if it does not exist, otherwise leave it untouched
open(os_path.join(main_path,param.sval('PARAMETER_CORRELATION_FILE')), 'a').close()
#
# create all iteration step paths, if they do not exist
for i in range(n_iter_steps):
    iter_path = os_path.join(main_path,param.sval('ITERATION_STEP_PATH')+'{0:03d}'.format(i+1))
    #
    if os_path.exists(iter_path):
        print "# NOTE, iteration step path '"+iter_path+"' already exists. Continuing, doing nothing."
    else:           
        create_iter_dir(iter_path,i+1)
print "###"
print ""
#------------------------------------------------------------
exit_now()
#------------------------------------------------------------
