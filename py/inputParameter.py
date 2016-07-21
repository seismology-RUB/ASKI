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
class inputParameter:
    """
    """
    def __init__(self,filename):
        self._key_val = dict( [
                ( line.split('=')[0].strip()  ,  line.split('=')[1].split('#')[0].strip() )
                for line in open(filename,'r')
                if line.strip()!=''            # ignore empty lines or lines containing spaces only
                if line.strip()[0:1]!='#'      # ignore comment lines (first non space character is '#')
                if '=' in line.split('#')[0]  ] ) # ignore lines which have no '=' in front of '#' (also works if not '#' in line)

    def hasKey(self,key):
        #return key in self._key_val
        return self._key_val.has_key(key)

    def keysNotPresent(self,keys):
        l = []
        for key in keys:
            #if not key in self._key_val:
            if not self._key_val.has_key(key):
                l.append(key)
        return l

    def prnt(self):
        for item in self._key_val.items():
            print "'"+"' = '".join(item)+"'"
            
    def getval(self,key):
        if self._key_val.has_key(key):
            return self._key_val[key]
        else:
            return None

    def ival(self,key):
        if self._key_val.has_key(key):
            try:
                return int(self._key_val[key])
            except:
                return None
        else:
            return None
        
    def ilist(self,key,n):
        if self._key_val.has_key(key):
            try:
                l = [int(i) for i in self._key_val[key].split()]
                if len(l) >= n:
                    return l[:n]
                else:
                    return None
            except:
                return None
        else:
            return None
        
    def fval(self,key):
        if self._key_val.has_key(key):
            try:
                return float(self._key_val[key])
            except:
                return None
        else:
            return None
        
    def flist(self,key,n):
        if self._key_val.has_key(key):
            try:
                l = [float(i) for i in self._key_val[key].split()]
                if len(l) >= n:
                    return l[:n]
                else:
                    return None
            except:
                return None
        else:
            return None
        
    def sval(self,key):
        if self._key_val.has_key(key):
            try:
                return self._key_val[key]
            except:
                return None
        else:
            return None

    def slist(self,key,n):
        if self._key_val.has_key(key):
            try:
                l = self._key_val[key].split()
                if len(l) >= n:
                    return l[:n]
                else:
                    return None
            except:
                return None
        else:
            return None
        
    def lval(self,key):
        if self._key_val.has_key(key):
            if self._key_val[key].lower() in ['.true.','true','t','1','1.']:
                return True
            elif self._key_val[key].lower() in ['.false.','false','f','0','0.']:
                return False
            else:
                return None
        else:
            return None
