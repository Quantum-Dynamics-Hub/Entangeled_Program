#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import generic_operations_ND


def compute_ETHD_pot_force(q, w, m, dim, approx, model):
#Returns the potential energy and force for the ETHD model, depending on the model and approximation chosen by the user.
#Params in:  q = list of trajectory positions
#            m = mass
#            approx = current approximation, ex) H = h0 or H = h0 + h1
#            model = the potential chosen by the user. Ex) 1 - cubic

    N = len(q)
    pot = [ [0.0]*dim, [0.0]*dim ]
    f = []

    if model == 1:
        pot, f = generic_operations_ND.harmonic_oscillator(q, w, m, dim)
    if model == 2:
        pot, f = generic_operations_ND.double_well_potential(q, dim)

    if approx == 2:
 
        # Computing H1
        # For a single dof, H1 = [ hbar^2 / 8.0*m*s_(dof)^2 ]
        # Where dof = alpha, s_alpha ^ 2 = (1/N) * SUM (q_alpha_i - q_alpha_avg)^2 
       
        # For each dof, we need a s^2 term:
        s2 = []
        # Each s^22 term needs the following components:
        sumq, sumqq, avg_q = [], [], []

        for j in range(dim):
            sumq.append( 0.0 )
            sumqq.append( 0.0 )
            avg_q.append( 0.0 )
            for i in range(N):
                # for N dof
                sumq[j] = sumq[j] + q[i][j] 
                sumqq[j] = sumqq[j] + q[i][j]*q[i][j]
                avg_q[j] = sumq[j]/float(N)

        # Now, sumq[N] is the sum of positions in the Nth dof.
        # We can now constrcut a s^2 term for each dof:
        
        h1 = []  
        m_avg = sum(m)/float(N)
        for j in range(dim):
            s2.append( (1.0/float(N))*sumqq[j] - (1.0/float(N)*1.0/float(N))*sumq[j]*sumq[j] )

            # The Perturbation term for each dof 
            pot[1][j] = pot[1][j] + N/(8.0*m_avg*s2[j])
        
        # Now, for each dof, let's calculate the force resulting from the perturbation term, H1:
        denom = []
        for j in range(dim):
            denom.append( 4.0*m_avg*s2[j]*s2[j] )

        for i in range(N):
            f.append( [0.0]*dim )

            for j in range(dim):
                f[i][j] = f[i][j] + (q[i][j] - avg_q[j])/denom[j]  

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
