#!/usr/bin/env python
#
#----------------------------------------------------------------------------
#   Copyright 2016 Wolfgang Friederich and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
"""
Python script that prints the content of file rules.mk, assuming to be called from
the ASKI installation directory (in which the src directory ./f90 is contained).
"""
import sys
import glob
import re
import os.path
#-----------------------------------------------------------
# define a list of modules for which "use module" statements should be ignored in .f90 files
ignore_modules_f90 = ['iso_fortran_env','mpi','iso_c_binding']
#-----------------------------------------------------------
def search_text_file(filename, pattern):
    """ Search each line of a text file for a pattern.

    Start searching at the beginning of the line.
    Returns the group(1) match.
    """
    matchlist = []
    f = open(filename,'r')
    for line in f:
        m = re.match(pattern,line)
        if m is not None:
            matchlist.append(m.group(1))
    f.close()
    return matchlist
#-----------------------------------------------------------
def build_makefile_entry(filename, matchlist, depext):
    """Build a makefile entry.

    Form: filename_without_path_and_extension.o: dependencies.
    matchlist: list of dependencies
    depext: extension to be appended to the dependencies
    """
    target = str.split(os.path.basename(filename),'.')[0]+'.o:'
    entry = target
    for el in set(matchlist):
         entry = entry+' '+el+depext
    return entry
#-----------------------------------------------------------
dirs = ['.','f90']   # list of directories including current one and those given on command line
#
#  search for use statements in f90 files
#
for folder in dirs:                                       # loop through folders
    for filename in glob.glob(folder+'/*.f90'):              # loop through Fortran 90 files
        matchlist = search_text_file(filename,'[\t ]*use (\w+)') # search for use ... patterns which have arbitrary leading number of space and tab
        # exclude f90 modules listed in ignore_modules_f90
        matchlist_excluded = [match for match in matchlist if not match in ignore_modules_f90]
        if len(matchlist_excluded)>0:
            entry = build_makefile_entry(filename,matchlist_excluded,'.o')
            print entry
