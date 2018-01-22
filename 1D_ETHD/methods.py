#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import generic_operations

def compute_RPMD_pot_force(q, w, m, approx, model):
#Returns the potential energy and force for the RPMD model, depending on the model and approximation chosen by the user.
#Params in:  q = list of trajectory positions
#            w = stiffness of spring
#            m = mass
#            approx = current approximation, ex) H = h0 or H = h0 + h1
#            model = the potential chosen by the user. Ex) 1 - double_well

    N = len(q)

    pot, f = [0.0, 0.0], [0.0]*N
    if model == 1:
        pot, f = generic_operations.cubic_potential(q, m)
    if model == 2:
        pot, f = generic_operations.double_well_potential(q)

    if approx == 2:       

        # For Potential
        for i in range(N):

            if i != N-1:
                pot[1] += (q[i] - q[i+1])*(q[i] - q[i+1])
            else:
                pot[1] += (q[i] - q[0])*(q[i] - q[0])

       # For Forces
        for i in range(N):

            if i == 0:
                f[i] = f[i] - m*w*w*(2.0*q[i] - q[N-1] - q[i+1])

            if i == N-1:
                f[i] = f[i] - m*w*w*(2.0*q[i] - q[i-1] - q[0])

            else:
                f[i] = f[i] - m*w*w*(2.0*q[i] - q[i-1] - q[i+1])

    pot[1] = pot[1] * 0.5*m*w*w

    return pot, f

def compute_ETHD_pot_force(q, m, approx, model):
#Returns the potential energy and force for the ETHD model, depending on the model and approximation chosen by the user.
#Params in:  q = list of trajectory positions
#            m = mass
#            approx = current approximation, ex) H = h0 or H = h0 + h1
#            model = the potential chosen by the user. Ex) 1 - cubic

    N = len(q)

    pot = [0.0, 0.0]; f = [0.0]*N
    if model == 1:
        pot, f = generic_operations.cubic_potential(q, m)
    if model == 2:
        pot, f = generic_operations.double_well_potential(q)

    if approx == 2:
 
        # Computing H1 = [ hbar^2 / 8.0*m*s^2 ]
       
        v, sumq, sumqq = 0.0, 0.0, 0.0
        for i in range(N):
            sumq += q[i]
            sumqq += q[i]*q[i] 

        sumq2 = sumq*sumq
        s2 = (sumqq/float(N)) - (sumq2/(float(N)*float(N)))
        s4 = s2*s2

        denom = 4.0*m*s4

        pot[1] = pot[1] + N/(8.0*m*s2)

        qavg = sumq / float(N)
        for i in range(N):
            f[i] = f[i] + (q[i] - qavg)/denom  

    return pot, f
