#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys                               
import os
import math

def harmonic_oscillator(q,m,w):

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
    # the form is: pot = [ [x_classic is here , y_classic is here] , [0.0, 0.0] ]
    # f will be a list of two component lists: [ [0,0],[1,1],[2,2],...,[N,N] ] 
    
    N = len(q)
    pot = [ [0.0,0.0], [0.0,0.0] ]
    f = []
    for i in xrange(N):

        # Potential: for the x dimension of each particle
        pot[0][0] = pot[0][0] + 0.5*m[i]*w*w*q[i][0]*q[i][0]

        # Potential: for the y dimension of each particle 
        pot[0][1] = pot[0][1] + 0.5*m[i]*w*w*q[i][1]*q[i][1]

        # Force: for the [x,y] dimension of each particle
        # Recall q[particle_index][dimension] -> q[0][1] = first particle, y-dimension
        f.append( [-m[i]*w*w*q[i][0],-m[i]*w*w*q[i][1]]  )

    #print "pot = ", pot
    #print "f = ", f
    #print "\n"
    #sys.exit(0)

    return pot, f

def compute_generic_pot_force(q,model):

    # 1D
    # pot will be [0.0, 0.0], where pot[0] holds the sum of classical contributions
    # f will be a list: [0,1,2,...,N] of forces, where each element of the list is the force for a particular particle

    # 2D
    # pot will be [ [0.0, 0.0] , [0.0, 0.0]  ] , where pot[0] holds the sum classical contributions
    # the form is: pot = [ [x_classic is here , y_classic is here] , [0.0, 0.0] ]
    # f will be a list of two component lists: [ [0,0],[1,1],[2,2],...,[N,N] ] 

    pot = [ [0.0,0.0], [0.0,0.0] ]
    f = []

    # Here we call the function that computes the generic potential and force
    # for a specific potential 
    if model == 1:
        pot, f = harmonic_oscillator(q,m,w)
    
    return sum(pot[0]), f                 
    
def kin_en(p,m):

    """
    Defines the kinetic energy for a Hamiltonian
    1D: T(p[0]) = 0.5*m*(p[0]*p[0])
    2D: T(p_x,p_y) = 0.5*(p_x*p_x + p_y*p_y)/m
  or T(p[i][0],p[i][1])
    Params in:
    p = momentum of particle
    m = list of masses 
    """    

    # The kinetic energy is strcutred like the potential energy for 2D, except we will not have the additional list for quantum forces
    # The kinetic energy is structured as follows: kin = [ [0.0,0.0] ]
    # In the form [ [x_kinetic_here , y_kinetic_here] ]

    N = len(p)
    kin = [ [0.0,0.0] ]

    for i in xrange(N):

        # Kinetic: for the x dimension of each particle
        kin[0][0] = kin[0][0] + 0.5*p[i][0]*p[i][0]/m[i]

        # Kinetic: for the y dimension of each particle
        kin[0][1] = kin[0][1] + 0.5*p[i][1]*p[i][1]/m[i]

    return sum(kin[0])

def propagate_p(p, f, dt):
# Returns the updated momenta
    N = len(p)
    for i in range(N):

        # Update x dimension of momentum for ith particle 
        p[i][0] = p[i][0] + f[i][0]*dt

        # Update y dimension of momentum for ith particle
        p[i][1] = p[i][1] + f[i][1]*dt

    return p

def propagate_q(q, p, m, dt):
# Returns the updated positions
    N = len(q)
    for i in range(N):
        
        # Update x dimension of position for ith particle
        q[i][0] = q[i][0] + p[i][0]*dt/m[i]

        # Update y dimension of position for ith particle
        q[i][1] = q[i][1] + p[i][1]*dt/m[i]

    return q


def main(p,q,m,w,Nsnaps,Nsteps,dt,model):

    N = len(q); t = 0
    
    T, Tx, Ty = 0.0, 0.0, 0.0
    V, Vx, Vy = 0.0, 0.0, 0.0
    Ex, Ey = 0.0, 0.0    

    os.system("mkdir energy")
    e = open("energy/energy.txt", "w")

    V, f = compute_generic_pot_force(q,model)   

    for i in range(Nsnaps):
    
        for j in range(Nsteps):

            t = t + dt

            for k in range(N):

                propagate_p(p, f, 0.5*dt)

                propagate_q(q, p, m, dt)
                V, f = compute_generic_pot_force(q,model)   

                propagate_p(p, f, 0.5*dt)

        T =  kin_en(p,m)
        E = T + V

        e.write( " %8.5f %8.5f %8.5f %8.5f\n" % (t, T/float(N), V/float(N), E/float(N)) )

# System Parameters
model = 1

# Now let's make two dimensions
q, p, m  = [], [], []
M = 10
mass = 2000.0
for i in range(M):
    q.append( [-0.1-0.01*i,-0.2+0.01*i] )
    p.append( [0.0,0.0] )
    m.append( mass )

#print "q = ", q
#print "p = ", p
#print "m = ", m
#sys.exit(0)

w = 0.1

# Simulation Parameters
Nsnaps = 10000
Nsteps = 1
dt = 0.1

main(p,q,m,w,Nsnaps,Nsteps,dt,model)
