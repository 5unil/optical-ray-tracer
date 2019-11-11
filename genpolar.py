import numpy as np

def rtpairs(R, N):
    for i in range(len(R)):
        r=R[i]
        n=N[i]
        for j in range(n):
            t = j*2*numpy.pi/n
            yield r,t