#----------------------------------------------------------------------------
#   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
try:
    import start as start
    cubit = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass


class ecart_invgrid():
    def __init__(self):
        self.filename_nodes = 'node_coordinates'
        self.filename_hex8 = 'cell_connectivity_hex8'
        self.filename_tet4 = 'cell_connectivity_tet4'
        #
        self.hex8 = []
        self.tet4 = []
        #
        print "number of volumes = "+str(cubit.get_volume_count())
        print "number of hexes in entire mesh: cubit.get_hex_count = "+str(cubit.get_hex_count())
        print "number of tets in entire mesh: cubit.get_tet_count =  "+str(cubit.get_tet_count())
        #
        for volume in cubit.get_entities("volume"):
            print "VOLUME "+str(volume)+":"
            print "  number of elements: "+str(cubit.get_volume_element_count(volume))
            hex8_list = cubit.get_volume_hexes(volume)
            self.hex8 += hex8_list
            print "  number of hexes:    "+str(len(hex8_list))
            tet4_list = cubit.get_volume_tets(volume)
            self.tet4 += tet4_list
            print "  number of tets:     "+str(len(tet4_list))
        #
        print "number of hexes found while iterating over all volumes: "+str(len(self.hex8))
        print "number of tets found while iterating over all volumes: "+str(len(self.tet4))


    def write(self,path=''):
        if len(path) != 0:
            if path[-1] != '/': path+='/'
        self.write_node(path)
        self.write_hex8(path)
        self.write_tet4(path)


    def write_node(self,path):
        f = open(path+self.filename_nodes,'w')
        nodes = cubit.parse_cubit_list('node','all')
        nnodes = len(nodes)
        print "writing node coordinates file '"+path+self.filename_nodes+"' containing "+str(nnodes)+" nodes"
        f.write('%10i\n'%nnodes)
        #
        for node in nodes:
            xyz = cubit.get_nodal_coordinates(node)
            f.write('%20f %20f %20f\n'%xyz)
        f.close()


    def write_hex8(self,path):
        f = open(path+self.filename_hex8,'w')
        nhex8 = len(self.hex8)
        print "writing hex8 file '"+path+self.filename_hex8+"' containing "+str(nhex8)+" hexahedral cells"
        f.write('%10i\n'%nhex8)
        #
        for hex8 in self.hex8:
            nodes = cubit.get_connectivity('Hex',hex8)
            f.write('%10i %10i %10i %10i %10i %10i %10i %10i\n'%nodes)
        f.close()


    def write_tet4(self,path):
        f = open(path+self.filename_tet4,'w')
        ntet4 = len(self.tet4)
        print "writing tet4 file '"+path+self.filename_tet4+"' containing "+str(ntet4)+" tetrahedral cells"
        f.write('%10i\n'%ntet4)
        #
        for tet4 in self.tet4:
            nodes = cubit.get_connectivity('Tet',tet4)
            f.write('%10i %10i %10i %10i \n'%nodes)
        f.close()


def export2ASKI(path_ASKI_export='./'):
    ASKI_invgrid=ecart_invgrid()
    ASKI_invgrid.write(path=path_ASKI_export)
    print 'finished export ecart inversion grid to ASKI'


if __name__ == '__main__':
    path='.'
    export2ASKI(path)
