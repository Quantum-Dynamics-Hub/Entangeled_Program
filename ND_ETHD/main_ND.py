#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import os
import math
import generic_operations_ND
import methods_ND

#"""
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
#"""

def main(qq, pq, w, q_grid, p_grid, Nsteps, Nsnaps, m, dim, dt, approx, model, ent_type, vel):

    q = qq
    p = pq

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
        V, f = methods_ND.compute_ETHD_pot_force(q, w, m, dim, approx, model)

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
        
                generic_operations_ND.propagate_p(p, f, dim, dt*0.5)

                generic_operations_ND.propagate_q(q, p, m, dim, dt)

                if ent_type == 1:
                    V, f = methods_ND.compute_ETHD_pot_force(q, w, m, dim, approx, model)

                generic_operations_ND.propagate_p(p, f, dim, dt*0.5)

                N = len(q)


        """
        s.write(" %8.5f  %8.5f" % (t, V[1]/float(N)) ),
        s.write( "\n" )
        """

        T = generic_operations_ND.compute_kin(p, m, dim)
        

        E = T + V

        ###==Printing Stuff==###
        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(dim*N), V/float(dim*N), E/float(dim*N)) )

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
##### Here we initilize the simulation #####

# Initilize lists of position, momentum, and mass
qq = []; pq = []; m = []
# Declare Number of trajectories 
Ntraj = 100 
dim = 1

rnd = Random()
sigma_q, sigma_p, q_mean, p_mean = [], [], [], []
for i in xrange(dim):

    q_mean.append( -1.1 )	                      
    p_mean.append( 0.0 )                 
    sigma_q.append( 0.04 )
    sigma_p.append( 0.0 )               

#"""
print "q_mean = ", q_mean
print "p_mean = ", p_mean
print "sigma_q = ", sigma_q
print "sigma_p = ", sigma_p
#sys.exit(0)
#"""

# Decalre system variables
mass = 2000
w = 0.1

# Initilize trajectories
for i in xrange(Ntraj):
    qq.append( [0.0]*dim )
    pq.append( [0.0]*dim )
    for j in range(dim):

        # If using Libra:
        qq[i][j] = q_mean[j] + sigma_q[j] * rnd.normal()
        pq[i][j] = p_mean[j] + sigma_p[j] * rnd.normal()

    # else
    #qq.append( [-1.2 + 0.005*i]*dim )
    #pq.append( [0.0]*dim )

    m.append( mass )

print "q = ", qq
print "p = ", pq
print "\n"
#sys.exit(0)

q_grid = []; p_grid = []
for i in range(-150,200):
    q_grid.append(0.01*i)

for i in range(-100,500):
    p_grid.append(0.1*i)                                

ent_type = 1      # 1 - ETHD, 2 - RPMD, 3 - QHD2
model = 2         # 1 - harmonic oscillator, 2 - double well 
vel_rescale = 0   # 0 - no, 1 - yes
approx = 2        # 1 - H = H0,  2 - H = H0 + H1
dt = 0.1
Nsnap = 500
Nstep = 1
#                          
main(qq, pq, w, q_grid, p_grid, Nstep, Nsnap, m, dim, dt, approx, model, ent_type, vel_rescale)  
