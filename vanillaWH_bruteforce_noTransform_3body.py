"""vanillaWH_bruteforce.py: Wisdom-Holman integrator solving the EOM with an ODE solver"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint
import integrator_tools

def drift(sim_jacobi, sim, h):
    """
    Advance (Keplerian-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        p0 (float): negative of Hamiltonian at time 0
        h (float): timestep

    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each updated particle in Jacobi coodinates
    """
    
    r1_x, r1_y, r1_z, v1_x, v1_y, v1_z = sim_jacobi[1,:6]
    r2_x, r2_y, r2_z, v2_x, v2_y, v2_z = sim_jacobi[2,:6]
    m1, m2, m3 = sim[:3,6]
    m1_jacobi, m2_jacobi, m3_jacobi = sim_jacobi[:3,6]
    
    mu1 = m1*m2 / (m1 + m2)
    mu2 = m3 * (m1 + m2) / (m1 + m2 + m3)
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    t = np.array([0, h])
    initial_vector1 = [r1_x, r1_y, r1_z, m2_jacobi*v1_x, m2_jacobi*v1_y, m2_jacobi*v1_z]
    initial_vector2 = [r2_x, r2_y, r2_z, m3_jacobi*v2_x, m3_jacobi*v2_y, m3_jacobi*v2_z]
    
    sol1 = odeint(drift_ODE1, initial_vector1, t, rtol=1e-13, atol=1e-13, args=(m1, m2, m3, mu1, mu2))
    sol2 = odeint(drift_ODE2, initial_vector2, t, rtol=1e-13, atol=1e-13,args=(m1, m2, m3, mu1, mu2))
    
    # Update position and velocity with h*r' and h*p' in Jacobi coordinates
    
    r1_x = sol1[-1,0]
    r1_y = sol1[-1,1]
    r1_z = sol1[-1,2]
    v1_x = sol1[-1,3]/m2_jacobi
    v1_y = sol1[-1,4]/m2_jacobi
    v1_z = sol1[-1,5]/m2_jacobi
    
    r2_x = sol2[-1,0]
    r2_y = sol2[-1,1]
    r2_z = sol2[-1,2]
    v2_x = sol2[-1,3]/m3_jacobi
    v2_y = sol2[-1,4]/m3_jacobi
    v2_z = sol2[-1,5]/m3_jacobi
    
    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1_z, v1_x, v1_y, v1_z])
    sim_jacobi[2,:6] = np.array([r2_x, r2_y, r2_z, v2_x, v2_y, v2_z])
    
    return sim_jacobi

def drift_ODE1(vector, t, r1, r2, m1, m2, m3, mu1, mu2):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    # Norm of Jacobi coordinates
    r1 = np.sqrt(rx**2 + ry**2 + rz**2)
    
    rx_eqn = px/mu1
    ry_eqn = py/mu1
    rz_eqn = pz/mu1
    
    px_eqn = -m1*m2 * r1**(-3) * rx
    py_eqn = -m1*m2 * r1**(-3) * ry
    pz_eqn = -m1*m2 * r1**(-3) * rz
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt

def drift_ODE2(vector, t, r1, r2, m1, m2, m3, mu1, mu2):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the second non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    # Norm of Jacobi coordinates
    r2 = np.sqrt(rx**2 + ry**2 + rz**2)
    
    rx_eqn = px/mu2
    ry_eqn = py/mu2
    rz_eqn = pz/mu2

    px_eqn = -m3*(m1+m2) * r2**(-3) * rx
    py_eqn = -m3*(m1+m2) * r2**(-3) * ry
    pz_eqn = -m3*(m1+m2) * r2**(-3) * rz
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
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
