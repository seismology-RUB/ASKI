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
import numpy as np
import matplotlib.pyplot as plt
from cmath import phase as cmath_phase
#
#############################################################
# better solve this by regular expression?!
def get_column_spectrum_file(filename,icolumn):
    return np.array([[complex(float(lsplit.split(',')[0]),float(lsplit.split(',')[1].split(')')[0])) 
      for lsplit in line.split('(')[1:]]for line in open(filename,'r')],dtype=np.complex64)[:,icolumn-1]
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
#############################################################
help_str = """
   plotSpectrum.py

USAGE: please give at least 2 arguments:
[1] spectrum file
[2] column to plot (if not given, will use first column)

EXAMPLES:
plotSpectrum.py ./data_060404_220503_TUR8_UP
plotSpectrum.py ./synthetics_060404_220503_TUR8 2
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
if len(arguments)<=0:
    print help_str
    sys_exit()
spectrum_file = arguments[0]
if len(arguments)>=2:
    try:
        i_spectrum_column = int(arguments[1])
    except:
        print "### ERROR: second argument '"+arguments[1]+"' could not be converted to integer"
        print help_str
        sys_exit()
else:
    i_spectrum_column = 1

# read in column of spectrum file
spectrum = get_column_spectrum_file(spectrum_file,i_spectrum_column)

figure_title = "'"+spectrum_file+"' , column "+str(i_spectrum_column)

fig = plt.figure()  #plt.figure(figsize=(9.5,9))
#fig.subplots_adjust(left=0.08,right=0.97,bottom=0.07,hspace=0.43)
axAMP = fig.add_subplot(121)
axAMP.set_title('amplitude', size = 13)
axAMP.set_xlabel('index of line in file')
#axAMP.set_ylabel('spectral unit')
axAMP.plot(abs(spectrum),linewidth=2)  #axAMP.plot(np.arange(0,NSTEP_orig)*dt_orig*1000,seis_array,linewidth=3)

axPHASE = fig.add_subplot(122)
axPHASE.set_title('phase', size = 13)
axPHASE.set_xlabel('index of line in file')
axPHASE.plot(np.angle(spectrum),linewidth=2)
axPHASE.set_ylim([-np.pi,np.pi])

fig.suptitle(figure_title, fontsize=14)

plt.show(block=True)
