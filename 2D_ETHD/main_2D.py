#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import os
import math
import generic_operations_2D
import methods_2D

#"""
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
#"""

def main(q, p, w, q_grid, p_grid, Nsteps, Nsnaps, m, dt, approx, model, ent_type, vel):

    N = len(q)
    t = 0.0; orig = float(N); barrier = 0.0; QG = len(q_grid); PG = len(p_grid) 
             
    """
    os.system("mkdir p_q_info")
    init_q = open("p_q_info/initial_q.txt","w")
    init_p = open("p_q_info/initial_p.txt","w")
    for i in range(N):
        init_q.write( " %8.5f" % (q[i]) )
        init_q.write("\n")
        init_p.write( " %8.5f" % (p[i]) )
        init_p.write("\n")

    w = generic_operations.compute_omega(q, m, approx, model)
    print "omega = ", w
    """

    if model == 1:
        barrier = 0.67
        q_max = 2.0
    if model == 2:
        barrier = 0.0
        q_max = 50.0

    """
    print "Calling RPMG force, see if V[1] Matches"
    print "V[0] = ", V[0]/float(N)
    print "V[1] = ", V[1]/float(N)
    """

    os.system("mkdir energy")
    e = open("energy/energy.txt", "w")

    """
    os.system("mkdir phase_space")
    os.system("mkdir pos_space")
    os.system("mkdir distribution_data")
    os.system("mkdir tunnel")
    os.system("mkdir perturbation")
    os.system("mkdir input")

    r = open("phase_space/phase_space.txt", "w")
    g = open("pos_space/pos_space.txt", "w")
    h = open("distribution_data/dist.txt", "w")
    hh = open("tunnel/tunnel.txt", "w")
    s = open("perturbation/perturbation.txt", "w")

    for i in range(N):
        r.write(" %8.5f  %8.5f" % (q[i], p[i]) ),
    r.write( "\n" )

    g.write( " %8.5f" % (t) ),
    for i in range(N):
        g.write( " %8.5f" % (q[i]) )
    g.write( "\n" )

    s.write(" %8.5f  %8.5f" % (t, V[1]/float(N)) ),
    s.write( "\n" )

    ### Printing Initial Distribution Information
    for j in range(QG):
        count = 0.00
        for i in range(N):
            if q[i] > q_grid[j-1] and q[i] < q_grid[j]:  
                count += 1.00
            prob = count/float(N)             
        h.write( str(q_grid[j]) + " " + str(prob) + "\n") 
    """

    ###==Printing Stuff==###
    #e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(N), V/float(N),E/float(N)) )

    # Initilize forces
    if ent_type == 1:
        V, f = methods_2D.compute_ETHD_pot_force(q, w, m, approx, model)

    # Propagate for Nsnaps
    for k in xrange(Nsnaps):

        """
        a = open("input/2D_phase_"+str(k)+".dat","w")
        for j in range(QG):
            for l in range(PG):
                count = 0.00
                for i in range(N):
                    if q[i] > q_grid[j-1] and q[i] < q_grid[j]:
                       if p[i] > p_grid[l-1] and p[i] < p_grid[l]:
                           count += 1.00
                prob = count/float(N)
                a.write(" %8.5f  %8.5f  %8.5f" % (q_grid[j], p_grid[l], prob) )
                a.write("\n")
            a.write("\n")
        """

        # Propagate for Nsteps
        for j in xrange(Nsteps):

            t = t + dt

            for k in xrange(N):
        
                generic_operations_2D.propagate_p(p, f, dt*0.5)

                generic_operations_2D.propagate_q(q, p, m, dt)

                if ent_type == 1:
                    V, f = methods_2D.compute_ETHD_pot_force(q, w, m, approx, model)

                generic_operations_2D.propagate_p(p, f, dt*0.5)

                N = len(q)


        """
        s.write(" %8.5f  %8.5f" % (t, V[1]/float(N)) ),
        s.write( "\n" )
        """

        T = generic_operations_2D.compute_kin(p, m)
        

        E = T + V

        ###==Printing Stuff==###
        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(2*N), V/float(2*N), E/float(2*N)) )

        """
        for i in range(N):
            r.write(" %8.5f  %8.5f" % (q[i], p[i]) ),
        r.write( "\n" )

        g.write( " %8.5f" % (t) ),
        for i in range(N):
            g.write( " %8.5f" % (q[i]) )
        g.write( "\n" )

        count = 0.0
        for i in range(N):
            if q[i] < barrier:
                count = count + 1.0
        count = count/orig
        hh.write( str(t) + " " + str(count) + " " + "\n")
        """
   
##########  End of Main Function  ##########        
############################################


rnd = Random()
# For positions
q1_mean = -1.1	                      
q2_mean = -1.1
sigma_q1 = 0.04
sigma_q2 = 0.04
# For momenta
p1_mean = 0.0                 
p2_mean = 0.0
sigma_p1 = 0.0
sigma_p2 = 0.0

q = []; p = []; m = []
M = 100 
mass = 2000
w = 0.1
for i in range(M):
    q.append( [q1_mean + sigma_q1 * rnd.normal(), q2_mean + sigma_q2 * rnd.normal()] )
    p.append( [p1_mean + sigma_p1 * rnd.normal(), p2_mean + sigma_p2 * rnd.normal()] )
    m.append( mass )
#    q.append( [-1.3+0.001*i , -1.3+0.001*i] )
#    p.append( [0.0,0.0] )

q_grid = []; p_grid = []
for i in range(-150,200):
    q_grid.append(0.01*i)

for i in range(-100,500):
    p_grid.append(0.1*i)                                

ent_type = 1      # 1 - ETHD, 2 - RPMD, 3 - QHD2
model = 2         # 1 - harmonic oscillator, 2 - double well 
vel_rescale = 0   # 0 - no, 1 - yes
approx = 2        # 1 - H = H0,  2 - H = H0 + H1
dt = 1.0
Nsnap = 100
Nstep = 1
#                          
main(q, p, w, q_grid, p_grid, Nstep, Nsnap, m, dt, approx, model, ent_type, vel_rescale)  
