"""
Python project aimed at building a simple 3-D optical ray tracer.
It can be used to model the behaviour of simple optical systems.
"""

import numpy as np

class Ray:
    def __init__(self, x=0, y=0, z=0, dx=0, dy=0, dz=0, p=None, k=None):
        self.__pAll = []
        self.__kAll = []
        Ray.append(self, x, y, z, dx, dy, dz, p, k)


    def currentPos(self):
        return self.__pAll[-1] ## Assuming last element added to list is current position


    def currentDir(self):
        return self.__kAll[-1] ## Assuming last element added to list is current direction


    def append(self, x=0, y=0, z=0, dx=0, dy=0, dz=0, p=None, k=None):
        if p==None:
            self.__p = np.array([x,y,z])
            self.__pAll.append(self.__p)
        else:
            if len(p) != 3:
                raise Exception("The position vector must have 3 components")
            else:
                self.__p = np.array(p)
                self.__pAll.append(self.__p)

        if k==None:
            self.__k = np.array([dx,dy,dz])
            self.__kAll.append(self.__k)
        else:
            if len(k) != 3:
                raise Exception("The direction vector must have 3 components")
            else:
                self.__k = np.array(k)
                self.__kAll.append(self.__k)


    def vertices(self):
        return "All positions: %r" % self.__pAll


    def __repr__(self):
        return "%s(p=%r, k=%r)" % ("Ray", self.__p, self.__k)
    

    def __str__(self):
        return "Position: %s, Direction: %s" % (self.__p, self.__k)
