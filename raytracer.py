"""
Python project aimed at building a simple 3-D optical ray tracer.
It can be used to model the behaviour of simple optical systems.
"""

import numpy as np

class Ray:
    """
    A ray class representing an optical ray using inputted positions and directions 
    """   
    def __init__(self, x=0, y=0, z=0, dx=0, dy=0, dz=0, p=None, k=None):
        self.__pAll = []
        self.__kAll = []
        self.append(x, y, z, dx, dy, dz, p, k)
    

    def currentPos(self):
        return self.__pAll[-1] # Assuming last element added to list is current position


    def currentDir(self):
        return self.__kAll[-1] # Assuming last element added to list is current direction


    def append(self, x=0, y=0, z=0, dx=0, dy=0, dz=0, p=None, k=None):
        if p==None:
            self.__p = np.array([x,y,z])
            self.__pAll.append(self.__p)
        elif len(p) != 3:
            raise Exception("The position vector must have 3 components")
        else:
            self.__p = np.array(p)
            self.__pAll.append(self.__p)

        if k==None:
            ktemp = np.array([dx,dy,dz])
            self.__k = ktemp/np.linalg.norm(ktemp) # Normalised k
            self.__kAll.append(self.__k)
        elif len(k) != 3:
            raise Exception("The direction vector must have 3 components")
        else:
            ktemp = np.array(k)
            self.__k = ktemp/np.linalg.norm(ktemp) # Normalised k
            self.__kAll.append(self.__k)


    def vertices(self):
        return "All positions: %r" % self.__pAll


    def __repr__(self):
        return "%s(p=%r, k=%r)" % ("Ray", self.__p, self.__k)

    
    def __str__(self):
        return "Position: %s, Normalised Direction: %s" % (self.__p, self.__k)


  
class OpticalElement:
    def propagate(self, ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()



class SphericalRefraction(OpticalElement):
    """
    z0 : Intersection of lens with z plane
    R : Aperture Radius
    curvature : 1/R
    n1 : Refractive index of first medium
    n2 : Refractive index of second medium
    """
    def __init__(self, z0, R, curv, n1, n2):
        self.__z0 = z0
        self.R = R # R is aperture radius not radius of curvature                    
        self.__curv = curv
        self.__n1 = n1
        self.__n2 = n2


    def intercept(self, ray):
        p = ray.currentPos()
        k = ray.currentDir()
        zmag = (self.__z0-p[2])/k[2] # magnitude of PQ vector found using z component analysis (Q = point on circumference)

        r = p-np.array([ 0, 0, self.__z0+(1/self.__curv)]) # r is the vector from centre of lens (not origin) to p.
        rsqr = (np.linalg.norm(r))**2
        rdotk = np.dot(r,k)  
        
        if self.__curv==0:
            q = np.array([p[0]+k[0]*zmag, p[1]+k[1]*zmag, self.__z0]) 
        elif (rdotk**2 < rsqr-(self.__curv**-2)):
            raise Exception("The ray does not intercept with the lens")
            return None
        elif self.__curv>0:
            l = -rdotk - ((rdotk)**2 - (rsqr-(self.__curv**-2))**0.5) # smaller valued solution to quadratic equation for l
            q = p+l*k
        else:
            l = -rdotk + ((rdotk)**2 - (rsqr-(self.__curv**-2))**0.5) # larger valued solution to quadratic equation for l
            q = p+l*k

        if abs((q[0]) > self.R) or abs((q[1]) > self.R):
            raise Exception("The ray does not intercept within the region of the lens")
            return None
        else:
            return q

        
    def refract(self, q, k1, n1, n2):
        if self.__curv==0:   
            N = [0,0,-1]
        else:
            r = q-np.array([ 0, 0, self.__z0+(1/self.__curv)]) # r is the vector from centre of lens (not origin) to q.
            N = r/np.linalg.norm(r) # N is directional unit vector from centre of curvature to point on circumference

        if (1-np.dot(N,k1)**2)**0.5 > (n2/n1):
            raise Exception("Total Internal Reflection is taking place")
            return None
        else:
            if (np.dot(N,k1) < 0):
                k2 = (n1/n2)*k1-N*((n1/n2)*np.dot(N,k1)+(1-((1-(np.dot(N,k1)**2))*(n1/n2)**2))**(0.5))
            else:
                k2 = (n1/n2)*k1-N*((n1/n2)*np.dot(N,k1)-(1-((1-(np.dot(N,k1)**2))*(n1/n2)**2))**(0.5))
            return k2


    def propagate(self, ray):
        p = self.intercept(ray) # returns position of intercept on lens
        k = self.refract(p, ray.currentDir(), self.__n1, self.__n2) # returns new direction after refraction
        ray.append(p,k)

    
    def __repr__(self):
        return "%s(z0=%g, R=%g, curv=%g, n1=%g, n2=%g)" % ("Spherical Refraction", self.__z0, self.R, self.__curv, self.__n1, self.__n2)



class OutputPlane(OpticalElement):
    """
    A class defined to propagate rays to their intersection point with the output plane
    """   
    def __init__(self, z1):
        self.__z1 = z1


    def intercept(self, ray):
        p = ray.currentPos()
        k = ray.currentDir()
        zmag = (self.__z1-p[2])/k[2] # magnitude of PQ vector found using z component analysis (Q = point on circumference)
        q = np.array([p[0]+k[0]*zmag, p[1]+k[1]*zmag, self.__z1]) # Considering exit side of lens is a plane surface
        return q


    def propagate(self,ray):
        p = self.intercept(ray)
        k = ray.currentDir()
        ray.append(p,k)




























