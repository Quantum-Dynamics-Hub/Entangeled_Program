#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import methods_ND

#"""
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
#"""

def harmonic_oscillator(q, w, m, dim):

    """
    Defines the potential for the Harmonic Oscillator
    1D: V(q[0]) = 0.5*m*w*w*(q[0]*q[0])
    2D: V(q_x,q_y) = 0.5*m*w*w(q_x*q_x + q_y*q_y)
    Params in:
    q = position of particle
    m = list of masses 
    w = potential specfifc constant 
    """    

    # Initilize list for potential energy and force
    # 1D
    # pot will be [0.0, 0.0], where pot[0] holds the sum of classical contributions
    # f will be a list: [0,1,2,...,N] of forces, where each element of the list is the force for a particular particle

    # 2D
    # pot will be [ [0.0, 0.0] , [0.0, 0.0]  ] , where pot[0] holds the sum classical contributions
    # the form is: pot = [ [1st_dof_classic is here , 2nd_dof_classic is here] , [0.0, 0.0] ]
    # f will be a list of two component lists: [ [0,0],[1,1],[2,2],...,[N,N] ] 
    
    N = len(q)
    pot = [ [0.0]*dim , [0.0]*dim ]
    f = []

    #print "pot = ", pot
    #print "f = ", f
    #print "\n"
    #sys.exit(0)

    for i in xrange(N):
        f.append( [0.0]*dim )    
        for j in xrange(dim):

            pot[0][j] = pot[0][j] +  0.5*m[i]*w*w*q[i][j]*q[i][j]

            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][j] = -m[i]*w*w*q[i][j] 

    #print "pot = ", pot
    #print "f = ", f
    #print "\n"
    #sys.exit(0)

    return pot, f

def double_well_potential(q, dim):
# Defines the double_well potential. 
# V(q) = A * ( 0.25*q^4 - 0.5*q^2 )
#Params in:  q = list of trajectory positions


    C = 1.0
    a, b = 0.25, 0.5

    N = len(q)
    pot = [ [0.0]*dim , [0.0]*dim ]
    f = []

    q2, q3, q4 = [], [], []
    for i in xrange(N):

        f.append( [0.0]*dim )
        q2.append( [0.0]*dim )
        q3.append( [0.0]*dim )
        q4.append( [0.0]*dim )    

        for j in xrange(dim):

            q2[i][j] = q[i][j]*q[i][j]
            q3[i][j] = q2[i][j]*q[i][j]
            q4[i][j] = q3[i][j]*q[i][j]

            pot[0][j] = pot[0][j] + C*( a*q4[i][j] - b*q2[i][j] )

            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][j] =  -C*( 4.0*a*q3[i][j] - 2.0*b*q[i][j] )

    print "pot = ", pot
    print "f = ", f
    print "\n"
    #sys.exit(0)

    return pot, f


def compute_kin(p,m, dim):

    """
    Defines the kinetic energy for a Hamiltonian
    1D: T(p[0]) = 0.5*(p[0]*p[0])/m
    2D: T(p[i][0],p[i][1]) = 0.5*(p[i][0]*p[i][0] + p[i][1]*p[i][1])/m

    Params in:
    p = momentum of particle
    m = list of masses 
    """    

    # The kinetic energy is structured like the potential energy for 2D, except we will not have the additional list for quantum forces
    # The kinetic energy is structured as follows: kin = [ [0.0,0.0] ]
    # In the form [ [1st_dof_kinetic_here , 2nd_dof_kinetic_here] ]

    N = len(p)
    kin = [ [0.0]*dim ]

    for i in xrange(N):
        for j in xrange(dim):

            # Kinetic: for the x dimension of each particle
            kin[0][j] = kin[0][j] + 0.5*p[i][j]*p[i][j]/m[i]

    return sum(kin[0])

def propagate_p(p, f, dim, dt):
# Returns the updated momenta

    N = len(p)
    for i in range(N):
        for j in xrange(dim):

            p[i][j] = p[i][j] + f[i][j]*dt

    return p

def propagate_q(q, p, m, dim, dt):
# Returns the updated positions

    N = len(q)
    for i in range(N):
        for j in xrange(dim):
        
            q[i][j] = q[i][j] + p[i][j]*dt/m[i]

    return q
