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
################################################################
#  This is the Makefile for program xgoertzel (for GNU Make)
################################################################
#
#-----------------------------------------------------------------------
#  set the compiler and compiler flags
#
COMPILER = gfortran
#
FFLAGS = -O3
#FFLAGS = -O3 -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-------------------------------------------------------------
#
default: all
all: compile run
#
#---------------------------------------------------------------
#
clean:
	-rm -f xgoertzel
#
#----------------------------------------------------------------
#
compile: clean
	$(COMPILER) $(FFLAGS) -o xgoertzel goertzel.f90
#
#----------------------------------------------------------------
#
run:
	@echo
	@echo +++ RUNNING EXECUTABLE ./xgoertzel NOW +++
	@echo
	./xgoertzel
#
#----------------------------------------------------------------
#
