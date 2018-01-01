#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import methods

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def double_well_potential(q):
# Defines the double_well potential. 
# V(q) = A * ( 0.25*q^4 - 0.5*q^2 )
#Params in:  q = list of trajectory positions
#            w = stiffness of spring
#            m = mass

    N = len(q)

    C = 1.0
    a, b = 0.25, 0.5
    pot = [0.0, 0.0]; f = [0.0]*N

    for i in range(N):
        q2 = q[i]*q[i]
        q3 = q[i]*q2
        q4 = q2*q2

        pot[0] = pot[0] + C*( a*q4 - b*q2 )   

        f[i] = -C*( 4.0*a*q3 - 2.0*b*q[i])

    return pot, f

def cubic_potential(q, m):
# V(q) = B*q^2 - B*q^3
#Params in:  q = list of trajectory positions
#            w = stiffness of spring
#            m = mass

    N = len(q)
  
    b, ww = 0.2981, 0.01
    A, B = b/3.0, 0.5*m*ww*ww

    pot, f = [0.0, 0.0], [0.0]*N

    for i in range(N):

        q2 = q[i]*q[i]
        q3 = q[i]*q2
        q4 = q2*q2

        pot[0] = pot[0] + ( (B*(q2)) - (A*q3) )

        f[i] = -( ( 2.0*B*(q[i])) - (3.0*A*q2 ) )     

    return pot, f


def compute_kin(p, m):
#Returns the kinetic energy, depending on the model and approximation chosen by the user.
#Params in:  p = list of trajectory momenta
#            m = mass

    N = len(p)
    sumpp = 0.0
    for i in range(N):
        sumpp += p[i]*p[i]
    return sumpp*0.5 / m   

def propagate_p(p, f, dt):
# Returns the updated momenta
    N = len(p)
    for i in range(N):
        p[i] = p[i] + f[i]*dt
    return p

def propagate_q(q, p, m, dt):
# Returns the updated positions
    N = len(q)
    for i in range(N):
        q[i] = q[i] + p[i]*dt/m
    return q

def compute_omega(q, m, approx, model):
# A temperature will be chosne such that the RPMD and EHTD perturbation terms are initially equal

    N = len(q)

    rpmd_sum, ethd_term = 0.0, 0.0
    for i in range(N):

        if i != N-1:
            rpmd_sum += (q[i] - q[i+1])*(q[i] - q[i+1])
        if i == N-1:
            rpmd_sum += (q[i] - q[0])*(q[i] - q[0])

    V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    ethd_term = V[1]
    Kb = 0.0000031668

    num = math.sqrt( 2.0*ethd_term )
    denom2 = Kb*N*math.sqrt( rpmd_sum*m )
    T = num/denom2

    # Compute omega
    omega = N*T*Kb

    # With omega now computed, we should check to see if the RPMD and ETHD perturbations are equal
    rpmd_term = 0.5*m*omega*omega*rpmd_sum
    print "Initial Value of RPMD perturbation = ", rpmd_term/float(N)
    print "Initial Value of ETHD perturbation = ", ethd_term/float(N)
    print "RPMD Temperature = ", T

    return omega

def traj_absorb(q, p, q_max):
    """
    q - [list of floats] trajectory positions
    p - [list of floats] trajectory momenta
    q_max - [float] position of trajectory absorber
    """

    N = len(q)

    # Find out which trajectories are to be removed
    removal_list = []
    flag = 0
    for i in range(N):
        if q[i] >= q_max:
            removal_list.append(i)    
            flag = 1

    # Update the lists with q and p, removing those that
    # corresponds to the trajectories crossed the absorber
    res_q = []
    res_p = []

    for i in xrange(N):
        if i in removal_list:
            pass
        else:
            res_q.append(q[i])
            res_p.append(p[i])

    return res_q, res_p

