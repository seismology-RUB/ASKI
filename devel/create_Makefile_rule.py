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
"""
Python script that creates a Makefile rule for an ASKI executable given by
sys.argv[1], also accounting for library requirements of its dependencies.
"""
from sys import argv as sys_argv
from sys import exit as sys_exit
from string import split,join
#
#############################################################
def add_to_list_if_not_present(long_list,add_list):
    for word in add_list:
        if word not in long_list:
            long_list.append(word)
    return long_list
#
#############################################################
def line_break(input_string):
    # break lines after max_char_line characters (in make line break syntax)
    max_char_line = 100
    input_words = split(input_string)
    broken_lines = []
    while len(input_words)>0 :
        len_input_words = [len(w) for w in input_words]
        #
        # n is that number of leading items in the current list input_words, 
        # for which the number of characters sum up to at most max_char_line 
        # (including space separation and appending two characters ' \')
        n = max([i+1 for i in range(len(len_input_words)) if sum(len_input_words[:i+1]) +i+2 <= max_char_line])
        #
        # append bare line content to list broken_lines without any line break characters or appending ' \'
        broken_lines.append(join(input_words[:n]))
        #
        # remove the respective words from list input_words
        input_words = input_words[n:]
    #
    # create output character string by joining the lines with (1) characters ' \' (indicating in GNU make a line 
    # break is following), (2) line break character '\n' and (3) tab character '\t' (indicating in GNU make that 
    # the previous line is continued)
    return join(broken_lines,' \\\n\t')
#
#############################################################
def compiler_of_target(target):
    if target in ['solveCglsKernelSystem']:
        return '$(MPICOMPILER)'
    else:
        return '$(COMPILER)'
#
#############################################################
def required_libs_for_depobs(obs):
    libs_list = []
    blas = '$(BLAS)'
    lapack = '$(LAPACK)'
    blacs = '$(BLACS)'
    scalapack = '$(SCALAPACK)'
    #
    # integrationWeights.o requires LAPACK (and BLAS?)
    if 'integrationWeights.o' in obs:
        add_to_list_if_not_present(libs_list,[blas,lapack])
    #
    # ecartInversionGrid.o requires BLAS
    if 'ecartInversionGrid.o' in obs:
        add_to_list_if_not_present(libs_list,[blas])
    #
    # kernelFocus.o requires BLAS and LAPACK
    if 'kernelFocus.o' in obs:
        add_to_list_if_not_present(libs_list,[blas,lapack])
    #
    # kernelLinearSystem.o requires BLAS
    if 'kernelLinearSystem.o' in obs:
        add_to_list_if_not_present(libs_list,[blas])
    #
    # solveCglsKernelSystem.o requires BLAS
    if 'solveCglsKernelSystem.o' in obs:
        add_to_list_if_not_present(libs_list,[blas])
    #
    # parallelLinearSystem.o requires BLACS and SCALAPACK
    if 'parallelLinearSystem.o' in obs:
        add_to_list_if_not_present(libs_list,[blacs,scalapack])
    #
    # serialLinearSystem.o requires LAPACK (and BLAS?)
    if 'serialLinearSystem.o' in obs:
        add_to_list_if_not_present(libs_list,[blas,lapack])
    #
    return join(libs_list)
#
#
#############################################################
# START MAIN SCRIPT HERE
#############################################################
#
# initialize variables
try:
    target = sys_argv[1]
except:
    print(__doc__)
    sys_exit()
depobs = []  # list of all objects on which target depends
rules_mk_lines = []  # lines of file rules.mk
#
# get started with dependencies of target and at the same time read and remember file 'rules.mk'
for line in open('rules.mk','r'):
    rules_mk_lines.append(line)
    if target in split(line)[0]:  # if 'target.o:' is the first word of line
        depobs = add_to_list_if_not_present(depobs,split(line)[1:])
#
# now add dependencies of each of the objects so far added to depobs
for obj in depobs:
    for line in rules_mk_lines:
        if obj in split(line)[0]:
            depobs = add_to_list_if_not_present(depobs,split(line)[1:])
#
# print (broken )line(s) of target dependencies
print line_break(target+': %: %.o '+join(depobs))
#
# print line to compile the target
compiler = compiler_of_target(target)
libs = required_libs_for_depobs([target+'.o']+depobs)
print join(['\t'+compiler+' -o $(bindir)/$@ $(obstring)',libs])
#
