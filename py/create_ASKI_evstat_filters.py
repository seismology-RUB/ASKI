#!/usr/bin/env python
#
#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.1.
#
#   ASKI version 1.1 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.1 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
from sys import argv as sys_argv
from sys import exit as sys_exit
from os import path as os_path
from os import mkdir as os_mkdir
from os import makedirs as os_makedirs
from inputParameter import inputParameter
from readEventStationFile import eventList,stationList
#
#############################################################
def exit_now():
    print "###"
    print "Good Bye"
    print "###"
    print ""
    sys_exit()
#############################################################
help_str = """
   create_ASKI_evstat_filters.py

USAGE: please give 1 argument:
[1] main parmeter file of the inversion

EXAMPLE:
create_ASKI_evstat_filters.py ./main_parfile_project1
"""
#
#
#############################################################
# START MAIN SCRIPT HERE
#############################################################
#
#
# creating station filters for all valid ASKI components: CX,CY,CZ,N,S,E,W,UP,DOWN
components = ['CX','CY','CZ','N','S','E','W','UP','DOWN']
#
arguments = sys_argv[1:]
#
# check if all required arguments are given
if len(arguments)!=1:
    print help_str
    sys_exit()
parfile = arguments[0]
if not os_path.exists(parfile):
    print "### ERROR: file '"+parfile+"' does not exist"
    print help_str
    sys_exit()
#
print ""
print "###"
print "running script create_ASKI_evstat_filters.py now, creating station and event filters all containing the real number 1. = 1. + 0.*i for all frequencies"
print "main parameter file is '"+parfile+"'"
print ""
#
param = inputParameter(parfile)
noKeys = param.keysNotPresent(['MEASURED_DATA_NUMBER_OF_FREQ','FILE_EVENT_LIST','FILE_STATION_LIST',
                                'PATH_EVENT_FILTER','PATH_STATION_FILTER'])
if len(noKeys) > 0:
    print "### ERROR : the following keywords are required in parameter file '"+parfile+"':"
    print "### "+',  '.join(noKeys)
    print ""
    exit_now()
#
nfreq = param.ival('MEASURED_DATA_NUMBER_OF_FREQ')
if nfreq is None:
    print "### ERROR : the value of 'MEASURED_DATA_NUMBER_OF_FREQ' in main parfile is not an integer"
    exit_now()
if nfreq < 0:
    print "### ERROR : the value of 'MEASURED_DATA_NUMBER_OF_FREQ' in main parfile must be positive (is "+str(nfreq)+")"
    exit_now()
#
# read event and stations lists
if os_path.exists(param.sval('FILE_EVENT_LIST')):
    evlist = eventList(param.sval('FILE_EVENT_LIST'),list_type='standard')
else:
    print "### ERROR : the event list file '"+param.sval('FILE_EVENT_LIST')+"' given in the main parfile does not exist"
    exit_now()
nev = evlist.nev
if os_path.exists(param.sval('FILE_STATION_LIST')):
    statlist = stationList(param.sval('FILE_STATION_LIST'),list_type='standard')
else:
    print "### ERROR : the station list file '"+param.sval('FILE_STATION_LIST')+"' given in the main parfile does not exist"
    exit_now()
nstat = statlist.nstat
#
if nev==0 and nstat==0:
    print "### no valid events in event list and no valid stations in station list, so nothing to do" 
    exit_now()
#
# create default event and station filters, containing all (1.0,0.0) for all frequencies
#
# filter file:
content_filter_file = nfreq*'  ( 1.0 , 0.0 )\n'
#
# create event filter files:
print "creating event filter files"
if not os_path.exists(param.sval('PATH_EVENT_FILTER')):
    print "  mkdir '"+param.sval('PATH_EVENT_FILTER')+"'"
    os_makedirs(param.sval('PATH_EVENT_FILTER'))
else:
    print "  path event filter '"+param.sval('PATH_EVENT_FILTER')+"' exists"
for iev in range(nev):
    filename = os_path.join(param.sval('PATH_EVENT_FILTER'),'filter_'+evlist.events[iev]['evid'])
    print "  creating filter file '"+filename+"'"
    open(filename,'w').write(content_filter_file)
print ""

print "creating station filter files"
if not os_path.exists(param.sval('PATH_STATION_FILTER')):
    print "  mkdir '"+param.sval('PATH_STATION_FILTER')+"'"
    os_makedirs(param.sval('PATH_STATION_FILTER'))
else:
    print "  path station filter '"+param.sval('PATH_STATION_FILTER')+"' exists"
for istat in range(nstat):
    for comp in components:
        filename = os_path.join(param.sval('PATH_STATION_FILTER'),'filter_'+statlist.stations[istat]['staname']+'_'+comp)
        print "  creating filter file '"+filename+"'"
        open(filename,'w').write(content_filter_file)
print ""

exit_now()
