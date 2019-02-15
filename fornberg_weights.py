def weights(z,x,m):
#===========================================================================
#
#   Calculates FD weights. The parameters are:
#   z   location where approximations are to be accurate,
#   x   vector with x-coordinates for grid points,
#   m   highest derivative that we want to find weights for
#   c   array size m+1,lentgh(x) containing (as output) in 
#       successive rows the weights for derivatives 0,1,...,m.
#
#   Adapted from "Finite-Difference Method" Scholarpedia Article
#   http://dx.doi.org/10.4249/scholarpedia.9685
#
#===========================================================================

    import numpy as np

    n=np.shape(x)[0]; c=np.zeros([m+1,n]); c1=1; c4=x[0]-z; c[0,0]=1
    for i in range(1,n):
        mn=np.min([i,m])+1; c2=1; c5=c4; c4=x[i]-z
        for j in range(i):
            c3=float(x[i]-x[j]); c2=float(c2*c3)
            if j==i-1:
                c[1:mn,i]=c1*(np.arange(1,mn)*c[:mn-1,i-1]-c5*c[1:mn,i-1])/c2
                c[0,i]=-c1*c5*c[0,i-1]/c2
            c[1:mn,j]=(c4*c[1:mn,j]-np.arange(1,mn)*c[:mn-1,j])/c3
            c[0,j]=c4*c[0,j]/c3
        c1=c2

    #only return the highest derivative
    #c_last = c[-1,:]

    return c

#Test the Python code on your example
fd_coeffs = weights(0,range(-2,3),6)
print('fd_coeffs = '); print(fd_coeffs)
