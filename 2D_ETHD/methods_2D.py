#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import generic_operations_2D


def compute_ETHD_pot_force(q, w, m, approx, model):
#Returns the potential energy and force for the ETHD model, depending on the model and approximation chosen by the user.
#Params in:  q = list of trajectory positions
#            m = mass
#            approx = current approximation, ex) H = h0 or H = h0 + h1
#            model = the potential chosen by the user. Ex) 1 - cubic

    N = len(q)
    pot = [ [0.0,0.0], [0.0,0.0] ]
    f = []

    if model == 1:
        pot, f = generic_operations_2D.harmonic_oscillator(q,m,w)
    if model == 2:
        pot, f = generic_operations_2D.double_well_potential(q)

    if approx == 2:
 
        # Computing H1
        # For a single dof, H1 = [ hbar^2 / 8.0*m*s_(dof)^2 ]
        # Where dof = alpha, s_alpha ^ 2 = (1/N) * SUM (q_alpha_i - q_alpha_avg)^2 
       
        sumq_1 = 0.0
        sumq_2 = 0.0

        sumqq_1 = 0.0
        sumqq_2 = 0.0

        m_avg = 2000.0
        for i in range(N):
            
            # 1st dof
            sumq_1 = sumq_1 + q[i][0]
            sumqq_1 = sumqq_1 + q[i][0]*q[i][0]

            # 2nd dof
            sumq_2 = sumq_2 + q[i][1]
            sumqq_2 = sumqq_2 + q[i][1]*q[i][1]
 
        # Method 1 (Should be equivalent to Method 2, simply a different way of expressing s^2 )
        #"""
        q_1_avg = sumq_1/float(N)
        q_2_avg = sumq_2/float(N)

        s2_1 = (1.0/float(N))*sumqq_1 - (1.0/(float(N)*float(N)))*sumq_1*sumq_1
        s2_2 = (1.0/float(N))*sumqq_2 - (1.0/(float(N)*float(N)))*sumq_2*sumq_2

        s4_1 = s2_1*s2_1
        s4_2 = s2_2*s2_2

        pot[1][0] = pot[1][0] + N/(8.0*m_avg*s2_1)
        pot[1][1] = pot[1][1] + N/(8.0*m_avg*s2_2)

        denom_1 = 4.0*m_avg*s4_1
        denom_2 = 4.0*m_avg*s4_2
        for i in range(N):

            f[i][0] = f[i][0] + (q[i][0] - q_1_avg)/denom_1  
            f[i][1] = f[i][1] + (q[i][1] - q_2_avg)/denom_2  
        #"""
        
        # Method 2 (Should be equivalent to Method 1, simply a different way of expressing s^2 )
        """
        q_1_avg = sumq_1/float(N)
        q_2_avg = sumq_2/float(N)

        s2_1 = 0.0
        s2_2 = 0.0
        for i in range(N):
        
            s2_1 = s2_1 + (1.0/float(N))*(q[i][0] - q_1_avg)*(q[i][0] - q_1_avg)
            s2_2 = s2_2 + (1.0/float(N))*(q[i][1] - q_2_avg)*(q[i][1] - q_2_avg)

        s4_1 = s2_1*s2_1
        s4_2 = s2_2*s2_2

        pot[1][0] = pot[1][0] + N/(8.0*m_avg*s2_1)
        pot[1][1] = pot[1][1] + N/(8.0*m_avg*s2_2)

        denom_1 = 4.0*m_avg*s4_1
        denom_2 = 4.0*m_avg*s4_2
        for i in range(N):

            f[i][0] = f[i][0] + (q[i][0] - q_1_avg)/denom_1  
            f[i][1] = f[i][1] + (q[i][1] - q_2_avg)/denom_2  
        """

    return ( sum(pot[0]) + sum(pot[1]) ) , f
