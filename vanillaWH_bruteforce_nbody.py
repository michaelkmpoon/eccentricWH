"""vanillaWH_bruteforce.py: Wisdom-Holman integrator solving the EOM with an ODE solver"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import solve_ivp
import integrator_tools

def drift(sim_jacobi, sim, h):
    """
    Advance (Keplerian-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        #p0 (float): negative of Hamiltonian at time 0
        h (float): timestep

    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each updated particle in Jacobi coodinates
    """
    
    # Apply the drift step for each non-central object
    for i in range(1, len(sim_jacobi)):
        
        rx, ry, rz, vx, vy, vz = sim_jacobi[i,:6]
        total_mass = np.sum(sim[:i+1,6])

        # Solve ODE for r' and v' in Jacobi coordinates
        t = np.array([0, h])
        initial_vector = [rx, ry, rz, vx, vy, vz]
        sol = solve_ivp(drift_ODE, t, initial_vector, method='RK45', t_eval = t, rtol=1e-13, atol=1e-13, args=(total_mass,))

        # Update position and velocity in Jacobi coordinates
        rx = sol.y[0,-1]
        ry = sol.y[1,-1]
        rz = sol.y[2,-1]
        vx = sol.y[3,-1]
        vy = sol.y[4,-1]
        vz = sol.y[5,-1]
        
        sim_jacobi[i,:6] = np.array([rx, ry, rz, vx, vy, vz])
    
    return sim_jacobi

def drift_ODE(t, vector, total_mass):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, vx, vy, vz = vector
    
    # separation between the centre-of-mass of objects interior to the orbiting object, 
    # and the orbiting object
    separation = np.sqrt(rx**2 + ry**2 + rz**2)
    
    rx_eqn = vx
    ry_eqn = vy
    rz_eqn = vz
    
    vx_eqn = -total_mass * separation**(-3) * rx
    vy_eqn = -total_mass * separation**(-3) * ry
    vz_eqn = -total_mass * separation**(-3) * rz
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn]
    
    return dvectordt

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
