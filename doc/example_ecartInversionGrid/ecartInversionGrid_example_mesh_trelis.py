#!python

#----------------------------------------------------------------------------
#   Copyright 2015 Florian Schumacher, Marcel Paffrath (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.0.
#
#   ASKI version 1.0 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.0 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------

# this python script was used with Trelis 15.1 in order to produce 
# the example ecart inversion grid as shown in the ASKI manual (ASKI version 1.0)

import cubit
import cubit2ASKIecartInversionGrid
reload(cubit2ASKIecartInversionGrid)

cubit.init([""])
cubit.cmd('reset')

################
# GEOMETRY
################
cubit.cmd('brick x 15 y 15 z 20')
cubit.cmd('move vol 1 x 0 y 0 z 10')
cubit.cmd('create Cylinder height 10 radius 2.5 ')
cubit.cmd('move vol 2 x 0 y 0 z 5')
cubit.cmd('subtract volume 2 from volume 1')
cubit.cmd('webcut vol all with plane zplane offset 15.0 noimprint nomerge')

################
# MESHING
################
cubit.cmd('vol 3 scheme Tetmesh')
cubit.cmd('surface 10 11 size 5')
cubit.cmd('surface 22 12 19 20 21 size 7')
cubit.cmd('mesh vol 3')
cubit.cmd('vol 1 scheme Tetmesh')
cubit.cmd('vol 1 size 5')
cubit.cmd('mesh vol 1')

################
# EXPORT
################
cubit.cmd('compress all')
cubit2ASKIecartInversionGrid.export2ASKI('EXPORT_ecartInversionGrid_manual')

