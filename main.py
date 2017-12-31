#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import os
import math
import generic_operations
import methods

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def rescale(q, p, w, m, q_max, approx, model, vel, ent_type):

    # Compute pre-absorption energy:
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    T_old = generic_operations.compute_kin(p, m)
    E_old = sum(V) + T_old

    q_new, p_new = generic_operations.traj_absorb(q, p, q_max)


    # Compute post-absorption energy:
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    T_new = generic_operations.compute_kin(p_new, m)
    E_new = sum(V) + T_new

    # Lets compensate for the energy change by rescaling momenta
    # instead of the old kinetic energy (subtract), we have a new one (add)
    # E_old = E_new - T_new + alp^2 * T_new 
    # So: E_old - E_new = T_new * ( alp^2 - 1)
    # alp^2 = 1  +  (E_old - E_new)/T_new 

    alp2 = 1.0 + (E_old - E_new)/T_new    # Ask

    alp = 1.0
    if alp2>0.0:
        alp = math.sqrt(alp2)

    if vel == 0:
        alp = 1.0

    N = len(p_new)
    for i in range(N):
        p_new[i] =  alp*p_new[i] 


    return q_new, p_new


def propagate_methods(q, p, w, m, f, dt, approx, model, ent_type):
    ############################################################
    # -------- Quantum Propagation: Entangeled Verlet -------- #
    ############################################################

    """

    N = len(q)
    for i in range(N):
        p[i] = p[i] + f[i]*dt*0.5
        q[i] = q[i] + p[i]*dt/m
        V, f = compute_pot_force(q, w, m, approx, model)
        p[i] = p[i] + f[i]*dt*0.5


    """

    generic_operations.propagate_p(p, f, dt*0.5)

    generic_operations.propagate_q(q, p, m, dt)


    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)


    generic_operations.propagate_p(p, f, dt*0.5)

    return q, p

#############################################################################
def main(q, p, w, q_grid, p_grid, q_max, Nsteps, Nsnaps, m, dt, approx, model, ent_type, vel):

    N = len(q)
    t = 0.0; orig = float(N); barrier = 0.0; QG = len(q_grid); PG = len(p_grid)


    # Initilize forces
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
        barrier = 0.67
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
        barrier = 0.0
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    V = sum(V)
    T = generic_operations.compute_kin(p, m)
    E = V + T

    os.system("mkdir energy")
    os.system("mkdir phase_space")
    os.system("mkdir pos_space")
    os.system("mkdir distribution_data")
    os.system("mkdir tunnel")

    e = open("energy/energy.txt", "w")
    r = open("phase_space/phase_space.txt", "w")
    g = open("pos_space/pos_space.txt", "w")
    h = open("distribution_data/dist.txt", "w")
    hh = open("tunnel/tunnel.txt", "w")
    
    ### Printing Initial Distribution Information
    for j in range(QG):
        count = 0.00
        for i in range(N):
            if q[i] > q_grid[j-1] and q[i] < q_grid[j]:  
                count += 1.00
            prob = count/float(N)             
        h.write( str(q_grid[j]) + " " + str(prob) + "\n") 

    ###==Printing Stuff==###
    e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(N), V/float(N),E/float(N)) )

    for i in range(N):
        r.write(" %8.5f  %8.5f" % (q[i], p[i]) ),
    r.write( "\n" )

    g.write( " %8.5f" % (t) ),
    for i in range(N):
        g.write( " %8.5f" % (q[i]) )
    g.write( "\n" )

    # Propagate for Nsnaps
    for k in xrange(Nsnaps):

        # Propagate for Nsteps
        for j in range(Nsteps):

            q, p = propagate_methods(q, p, w, m, f, dt, approx, model, ent_type)

            t = t + dt

            q, p = rescale(q, p, w, m, q_max, approx, model, vel, ent_type)
            N = len(q)


            if ent_type == 1:
                V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
            if ent_type == 2:
                V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
            if ent_type == 3:
                V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

        
        V = sum(V)
        T = generic_operations.compute_kin(p, m)
        E = V + T

        ###==Printing Stuff==###
        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(N), V/float(N),E/float(N)) )

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


##########  End of Main Function  ##########        
############################################

rnd = Random()
sigma_q = 0.04
sigma_p = 0.0               
q_mean = -1.1	                      
p_mean = 0.0                 
M = 5000                           
q = []; p = []
for i in range(M):
    q.append(q_mean + sigma_q * rnd.normal() )
    p.append(p_mean + sigma_p * rnd.normal() )
#    q.append(-0.45 + 0.005*i)
#    p.append(0.0)

q_max = 50.0

q_grid = []; p_grid = []
for i in range(-200,500):
    q_grid.append(0.01*i)

for i in range(-100,500):
    p_grid.append(0.1*i)                                
#                                        Nstep     Nsnap       mass        dt       approx       model        Entanglement_type       vel_rescale
main(q, p, 1.0, q_grid, p_grid, q_max,    100,     4200,      2000.0,      0.1,       2,          2,                1,                    0)                                                      
