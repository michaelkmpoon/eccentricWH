"""vanillaWH.py: Minimal Wisdom-Holman integrator following WHFast"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
import integrator_tools

def drift(sim_jacobi, sim, h):
    """
    Advance (Keplerian part) system by timestep h
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        h (float): timestep

    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each updated particle in Jacobi coodinates  
    """
    
    # Apply the drift step for each non-central object
    for i in range(1, len(sim_jacobi)):
    
        rx, ry, rz, vx, vy, vz = sim_jacobi[i, :6]

        total_mass = np.sum(sim[:i+1,6])

        r0 = np.sqrt(rx**2 + ry**2 + rz**2)
        r0i = 1./r0
        v2 = vx**2 + vy**2 + vz**2
        beta = 2*total_mass*r0i - v2
        eta0 = rx*vx + ry*vy + rz*vz
        zeta0 = total_mass - beta*r0

        # Solve Kepler's equation for X
        dtr0i = h * r0i # first order guess
        X_guess = dtr0i * (1 - dtr0i*eta0*0.5*r0i) # second order guess (following REBOUND)

        X = optimize.newton(func=Keplers_eqn, x0=X_guess, args=(r0, beta, eta0, zeta0, h))

        G1 = calculate_G(1, beta, X) 
        G2 = calculate_G(2, beta, X)
        G3 = calculate_G(3, beta, X)

        ri = 1./(r0 + eta0*G1 + zeta0*G2)
        f = 1 - total_mass*G2*r0i
        g = h - total_mass*G3
        fd = -total_mass*G1*r0i*ri   # time derivative of f
        gd = 1 - total_mass*G2*ri    # time derivative of g

        # solve for new position and velocity using f and g functions
        sim_jacobi[i, 0] = f*rx + g*vx
        sim_jacobi[i, 1] = f*ry + g*vy
        sim_jacobi[i, 2] = f*rz + g*vz
        sim_jacobi[i, 3] = fd*rx + gd*vx
        sim_jacobi[i, 4] = fd*ry + gd*vy
        sim_jacobi[i, 5] = fd*rz + gd*vz

    return sim_jacobi

def Keplers_eqn(X, r0, beta, eta0, zeta0, h):
    """
    Kepler's equation as described in eqn. (11) of Mikkola & Innanen 1999
    """
    term1 = r0 * X
    term2 = eta0 * calculate_G(2, beta, X)
    term3 = zeta0 * calculate_G(3, beta, X)
    return term1 + term2 + term3 - h

def calculate_G(n, beta, X):
    """
    G-functions as described in eqn. (9) of Mikkola & Innanen 1999
    """
    return X**n * calculate_c(n, beta*X**2)

def c(n, z, j):
    """
    Helper function to calculate_c
    """
    return (-z)**j / (math.factorial(n + 2*j))
    
def calculate_c(n, z, tolerance=1e-14, max_j=30):
    """
    c-functions as described in eqn. (7) of Mikkola & Innanen 1999
    
    Parameters:
        n (int): input for c-function
        z (float): input for c-function
        tolerance (float): stop expansion when adding the jth term 
            changes the relative c by less than tolerance. 
            Default is 1e-14.
        j (int): max number of terms in c-functione expansion
            Default is 30.
    
    Returns:
        (float): c-function value at given n, z.
    """
    
    j = 2 # number of terms in series expansion, after 0 and 1
    
    current_c = c(n, z, 0)
    next_c = current_c + c(n, z, 1)
    
    while (abs((current_c-next_c)/current_c) > tolerance) and j < max_j:
        current_c = next_c
        next_c += c(n, z, j)
        j += 1
        
    if j >= max_j:
        print('warning: c-function not converged')
        
    return next_c

def kick(sim_jacobi, sim, h):
    """
    Advance (Interaction part) system by timestep h
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        h (float): timestep

    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each updated particle in Jacobi coodinates  
    """

    # (1/4) Calculate x, y, z accelerations from first part of the Interaction Hamiltonian (in Jacobi coords)

    acceleration1j = np.zeros((len(sim_jacobi), 3))
    
    Mi = 0.
    Mi += sim[0, 6] + sim[1, 6]

    for i in range(2, len(sim_jacobi)):

        mi = sim[i, 6]
        Mi += mi

        # x, y, z in jacobi coordinates for the ith object 
        rij_x = sim_jacobi[i, 0]
        rij_y = sim_jacobi[i, 1]
        rij_z = sim_jacobi[i, 2]
        rij3 = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)**3

        acceleration1j[i, 0] = Mi * (rij_x / rij3)
        acceleration1j[i, 1] = Mi * (rij_y / rij3) 
        acceleration1j[i, 2] = Mi * (rij_z / rij3) 

    # (2/4) Calculate k-th x, y, z accelerations from second part of the Interaction Hamiltonian (in inertial coords)
    # 
    # a_{k,inertial,x} = Sum (from k=0 to N-1, j=k+1 to N-1) of mk * mj * (rj_x - rk_x) / norm (rk - rj)^3
    #                  + Sum (from i=0 to N-1, k=i+1 to N-1) of mi * mk * (rk_x - rj_x) / norm (ri - rk)^3

    acceleration2 = np.zeros((len(sim_jacobi), 3))
    
    for k in range(len(sim_jacobi)):
        
        for j in range(k+1, len(sim_jacobi)):
            
            if j != 1:
                
                mj = sim[j, 6]

                # here, j is the jth term, not Jacobi coordinate
                rk_x, rj_x = sim[k, 0], sim[j, 0] 
                rk_y, rj_y = sim[k, 1], sim[j, 1] 
                rk_z, rj_z = sim[k, 2], sim[j, 2] 
                dx = rk_x - rj_x
                dy = rk_y - rj_y
                dz = rk_z - rj_z
                rk_minus_rj3 = np.sqrt(dx**2 + dy**2 + dz**2)**3

                acceleration2[k, 0] += mj * (rj_x-rk_x) / rk_minus_rj3
                acceleration2[k, 1] += mj * (rj_y-rk_y) / rk_minus_rj3
                acceleration2[k, 2] += mj * (rj_z-rk_z) / rk_minus_rj3
            
    for i in range(len(sim_jacobi)):
        
        for k in range(i+1, len(sim_jacobi)):
            
            if k != 1:

                mi = sim[i, 6]

                ri_x, rk_x = sim[i, 0], sim[k, 0] 
                ri_y, rk_y = sim[i, 1], sim[k, 1] 
                ri_z, rk_z = sim[i, 2], sim[k, 2] 
                dx = ri_x - rk_x
                dy = ri_y - rk_y
                dz = ri_z - rk_z
                ri_minus_rk3 = np.sqrt(dx**2 + dy**2 + dz**2)**3

                acceleration2[k, 0] -= mi * (rk_x-ri_x) / ri_minus_rk3
                acceleration2[k, 1] -= mi * (rk_y-ri_y) / ri_minus_rk3
                acceleration2[k, 2] -= mi * (rk_z-ri_z) / ri_minus_rk3
            
    # (3/4) Convert accelerations from second part of the Interaction Hamiltonian to Jacobi coords
    
    acceleration2j = integrator_tools.inertial_to_jacobi_acc(sim, acceleration2)
    
    # (4/4) kick Jacobi velocities by h * (accel_part1 + accel_part2)
    
    for i in range(len(sim_jacobi)):
        
        sim_jacobi[i, 3] += h * (acceleration1j[i, 0] + acceleration2j[i, 0]) 
        sim_jacobi[i, 4] += h * (acceleration1j[i, 1] + acceleration2j[i, 1]) 
        sim_jacobi[i, 5] += h * (acceleration1j[i, 2] + acceleration2j[i, 2]) 

    return sim_jacobi
    
