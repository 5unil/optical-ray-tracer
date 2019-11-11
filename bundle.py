import numpy
import raytracer
import matplotlib.pylab as pyl

def rtpairs(R, N):
    for i in range(len(R)):
        r=R[i]
        n=N[i]
        for j in range(n):
            t = j*2*numpy.pi/n
            yield r,t
            
def bundlerays(n, rmax, m):
    bund = []
    R=numpy.arange(0,rmax,rmax/n)
    N=numpy.arange(1, n*m, m) 
    for r,t in rtpairs(R, N):
        myRay = raytracer.Ray(r * numpy.cos(t), r * numpy.sin(t),-200.,0,0,1.)        
        bund.append(myRay)
    return bund

def work(n, rmax, m):
    q=0
    u=0
    for ray in bundlerays(n, rmax, m):
        sphere = raytracer.SphericalRefraction(2.,0.03,1.,1.5,20.)
        H=raytracer.OutputPlane(150)
        sphere.propagate(ray)
        H.propagate(ray)
        k = ray.vertices()
        tmp = zip(*k)
        pyl.plot(tmp[2],tmp[0])
        #pyl.plot(tmp[1][-1],tmp[0][-1],'bo')
        #pyl.plot(tmp[1][0],tmp[0][0],'ro')
        q += tmp[0][-1]**2 + tmp[1][-1]**2
        u +=1 
    #print "RMS Deviation", (q/u)**(0.5)
        
        

#runloopg(10,0.1,6)
