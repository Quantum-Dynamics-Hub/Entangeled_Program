#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import os
import math

"""
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
"""


def cubic_potential(q, s, m):
# V(q) = B*q^2 - B*q^3
# Params in: q = mean position
#            s = wavepacket width
#            f = force
#            m = mass
 
    A = 0.33
    B = 0.5

    # V(q)
    pot =  B*q*q - A*q*q*q 

    # 0.5*s*s*V''(q)
    pot2 = 0.5*s*s*(2.0*B - 6.0*A*q) 

    term = 1.0 / (8.0*m*s*s)

    f_q =  3.0*A*(q*q + s*s) - 2.0*B*q       
    f_s = (1.0 / (4.0*m*s*s*s)) - s*(2.0*B - 6.0*A*q)

    V = pot + pot2 + term

    return V, f_q, f_s

def double_well_potential(q, s, m):
# V(q) =      C * (0.25*q^4 - 0.5*q^2)
# V'(q) =     C * (q^3 - q)
# V''(q) =    C * (3.0*q^2 - 1.0)
# V'''(q) =   C * (6.0*q)
# V''''(q) =  C * (6.0)
# Params in: q = mean position
#            s = wavepacket width
#            f = force
#            m = mass

    C = 1.0

    # V(q)
    pot =  C*(0.25*q*q*q*q - 0.5*q*q) 

    # 0.5*s*s*V''(q)
    pot2 = 0.5*s*s*C*(3.0*q*q - 1.0) 

    # s^4 * (1/8) * V''''(q) 
    # V''''(q) = 6*C
    pot4 = s*s*s*s*0.75*C 

    term = 1.0 / (8.0*m*s*s)

    f_q =  -C*(q*q*q - q) - 0.5*s*s*6.0*C*q     
    f_s = (1.0 / (4.0*m*s*s*s)) - s*C*(3.0*q*q - 1.0) - 3.0*s*s*s*C

    V = pot + pot2 + pot4 + term

    return V, f_q, f_s


def compute_pot_force(q, s, f, f_s, m, model):
#Returns the potential energy and force, depending on the model and approximation chosen by the user.
# Params in: q = mean position
#            s = wavepacket width
#            f_p and f_ps = force
#            m = mass
#            model = the potential chosen by the user. Ex) 1 - cubic

    if model == 1:
        V, f, f_s = cubic_potential(q, s, m)       
    if model == 2:
        V, f, f_s = double_well_potential(q, s, m)       

    return V, f, f_s

def compute_kin(p, ps, m):
#Returns the kinetic energy, depending on the model and approximation chosen by the user.
#Params in:  p = mean momenta
#            ps = coupled position + momenta
#            m = mass

    term1 = p*p*0.5/m
    term2 = ps*ps*0.5/m 
    return term1 + term2

# Here we will define the propagation of the Hamiltonian #
def propagate_H(q, p, s, ps, f, f_s, m, dt, model):

    p = p + f*dt*0.5
    ps = ps + f_s*dt*0.5

    s = s + ps*dt/m
    q = q + p*dt/m

    V, f, f_s = compute_pot_force(q, s, f, f_s, m, model)       

    ps = ps + f_s*dt*0.5
    p = p + f*dt*0.5

    return q, p, s, ps 


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
def main(qq, pq, s, ps, m, Nsteps, Nsnaps, dt, model):

    t = 0.0; 

    f = 0.0
    f_s = 0.0
   
    q = qq
    p = pq

    os.system("mkdir energy")
    os.system("mkdir phase_space")
    os.system("mkdir pos_space")
    os.system("mkdir distribution_data")
    os.system("mkdir tunnel")
    os.system("mkdir s_space")

    e = open("energy/energy.txt", "w")
    r = open("phase_space/phase_space.txt", "w")
    b = open("s_space/s_space.txt", "w")
    g = open("pos_space/pos_space.txt", "w")
    h = open("tunnel/tunnel.txt", "w")


    # Initilize forces
    V, f, f_s = compute_pot_force(q, s, f, f_s, m, model)       

    for i in range(Nsnaps):

        for j in range(Nsteps):

            q, p, s, ps = propagate_H(q, p, s, ps, f, f_s, m, dt, model)                                                                                                    
            t = t + dt 

            V, f, f_s = compute_pot_force(q, s, f, f_s, m, model)       

       

        T = compute_kin(p, ps, m)
        E = T + V

        # 12/14/17 - We are trying to use the error function:
        prob =  0.5 * ( math.erf( (3.0-q)/(math.sqrt(2)*s)) - math.erf( (0.0-q)/(math.sqrt(2)*s) ) ) 
        
        ###==Printing Stuff==###
        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T, V, E) )

        r.write(" %8.5f  %8.5f" % (q, p) ),
        r.write( "\n" )

        b.write(" %8.5f  %8.5f" % (t, s) ),
        b.write( "\n" )

        g.write( " %8.5f  %8.5f" % (t, q) ),
        g.write( "\n" )

        h.write(" %8.5f  %8.5f" % (t, prob) ),
        h.write( "\n" )

##########  End of Main Function  ##########        
############################################

qq = -1.1
pq =  0.0
s = 0.04
ps = 0.0

#                     mass     Nsteps   Nsnaps      dt    model  
main(qq, pq, s, ps,   2000.0,   1,      100000,    0.001,    2) 
