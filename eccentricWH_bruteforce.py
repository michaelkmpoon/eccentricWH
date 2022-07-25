"""eccentricWH_bruteforce.py: Wisdom-Holman integrator (using ODE solver) for highly eccentric orbits following Mikkola 1997"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint
import vanillaWH_tools

def drift(sim_jacobi, sim, A1, A2, p0, h):
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
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(sim_jacobi[1,0]**2 + sim_jacobi[1,1]**2 + sim_jacobi[1,2]**2)
    r2 = np.sqrt(sim_jacobi[2,0]**2 + sim_jacobi[2,1]**2 + sim_jacobi[2,2]**2)
    
    mu1 = m1*m2 / (m1 + m2)
    mu2 = m3 * (m1 + m2) / (m1 + m2 + m3)
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    t = np.array([0, h])
    initial_vector1 = [r1_x, r1_y, r1_z, m2*v1_x, m2*v1_y, m2*v1_z]
    initial_vector2 = [r2_x, r2_y, r2_z, m3*v2_x, m3*v2_y, m3*v2_z]
    
    sol1 = odeint(drift_ODE1, initial_vector1, t, args=(r1, r2, A1, A2, m1, m2, m3, mu1, mu2, p0))
    sol2 = odeint(drift_ODE2, initial_vector2, t, args=(r1, r2, A1, A2, m1, m2, m3, mu1, mu2, p0))
    
    # Update position and velocity with h*r' and h*p' in Jacobi coordinates
    
    r1_x = sol1[-1,0]
    r1_y = sol1[-1,1]
    r1_z = sol1[-1,2]
    v1_x = sol1[-1,3]/m2
    v1_y = sol1[-1,4]/m2
    v1_z = sol1[-1,5]/m2
    
    r2_x = sol2[-1,0]
    r2_y = sol2[-1,1]
    r2_z = sol2[-1,2]
    v2_x = sol2[-1,3]/m3
    v2_y = sol2[-1,4]/m3
    v2_z = sol2[-1,5]/m3
    
    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1_z, v1_x, v1_y, v1_z])
    sim_jacobi[2,:6] = np.array([r2_x, r2_y, r2_z, v2_x, v2_y, v2_z])
    
    return sim_jacobi

def drift_ODE1(vector, t, r1, r2, A1, A2, m1, m2, m3, mu1, mu2, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    rx_eqn = time_transform(r1, r2, A1, A2)*px/mu1
    ry_eqn = time_transform(r1, r2, A1, A2)*py/mu1
    rz_eqn = time_transform(r1, r2, A1, A2)*pz/mu1
    
    p_squared = (px**2 + py**2 + pz**2)
    Kepler_term1 = p_squared/(2*mu1) - m1*m2/r1
    Kepler_term2 = p_squared/(2*mu2) - m3*(m1+m2) / r2
    Kepler_terms = Kepler_term1 + Kepler_term2 + p0
    
    dKepler_terms_wrt_r1x = m1*m2 * r1**(-3) * rx
    dKepler_terms_wrt_r1y = m1*m2 * r1**(-3) * ry
    dKepler_terms_wrt_r1z = m1*m2 * r1**(-3) * rz
    
    term1_x = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r1x
    term1_y = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r1y
    term1_z = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r1z
    
    term2_x = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * rx * r1**(-1.5)  
    term2_y = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * ry * r1**(-1.5)  
    term2_z = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * rz * r1**(-1.5)    
    px_eqn = - (term1_x + term2_x)
    py_eqn = - (term1_y + term2_y)
    pz_eqn = - (term1_z + term2_z)
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt

def drift_ODE2(vector, t, r1, r2, A1, A2, m1, m2, m3, mu1, mu2, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the second non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    rx_eqn = time_transform(r1, r2, A1, A2)*px/mu2
    ry_eqn = time_transform(r1, r2, A1, A2)*py/mu2
    rz_eqn = time_transform(r1, r2, A1, A2)*pz/mu2
    
    p_squared = (px**2 + py**2 + pz**2)
    Kepler_term1 = p_squared/(2*mu1) - m1*m2/r1
    Kepler_term2 = p_squared/(2*mu2) - m3*(m1+m2) / r2
    Kepler_terms = Kepler_term1 + Kepler_term2 + p0
    
    dKepler_terms_wrt_r2x = m3*(m1+m2) * r2**(-3) * rx
    dKepler_terms_wrt_r2y = m3*(m1+m2) * r2**(-3) * ry
    dKepler_terms_wrt_r2z = m3*(m1+m2) * r2**(-3) * rz
    
    term1_x = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r2x
    term1_y = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r2y
    term1_z = time_transform(r1, r2, A1, A2) * dKepler_terms_wrt_r2z
    
    term2_x = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * rx * r2**(-1.5)  
    term2_y = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * ry * r2**(-1.5)  
    term2_z = Kepler_terms * (-time_transform(r1, r2, A1, A2)**2) * rz * r2**(-1.5)    
    
    px_eqn = - (term1_x + term2_x)
    py_eqn = - (term1_y + term2_y)
    pz_eqn = - (term1_z + term2_z)
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt

def time_transform(r1, r2, A1, A2):
    return 1. / (A1/r1 + A2/r2) 

def time_transform_derivative1(r1, r2, A1, A2):
    """
    Derivative of time transform w.r.t. r1
    """
    return A1 * (A1/r1 + A2/r2)**(-2) * r1**(-2) 

def time_transform_derivative2(r1, r2, A1, A2):
    """
    Derivative of time transform w.r.t. r2
    """
    return A2 * (A1/r1 + A2/r2)**(-2) * r2**(-2) 

def kick(sim_jacobi, sim, A1, A2, h):
    """
    Advance (Interaction-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        p0 (float): negative of Hamiltonian at time 0
        h (float): timestep

    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each updated particle in Jacobi coodinates  
    """
    
    m1, m2, m3 = sim[:3,6]
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(sim_jacobi[1,0]**2 + sim_jacobi[1,1]**2 + sim_jacobi[1,2]**2)
    r2 = np.sqrt(sim_jacobi[2,0]**2 + sim_jacobi[2,1]**2 + sim_jacobi[2,2]**2)
    
    # (1/4) Calculate x, y, z accelerations from first part of the 
    #       Interaction-like Hamiltonian (in Jacobi coords)
    
    acceleration1j = np.zeros((len(sim_jacobi), 3))
    
    # Only the 2nd term for m3
    r2_x = sim_jacobi[2, 0]
    r2_y = sim_jacobi[2, 1]
    r2_z = sim_jacobi[2, 2]
       
    term1_x = time_transform(r1, r2, A1, A2) * (- (m1+m2)) * r2_x * r2**(-1.5) 
    term1_y = time_transform(r1, r2, A1, A2) * (- (m1+m2)) * r2_y * r2**(-1.5) 
    term1_z = time_transform(r1, r2, A1, A2) * (- (m1+m2)) * r2_z * r2**(-1.5) 
    term2_x = (m1+m2)/r2 * (-time_transform(r1, r2, A1, A2)**2) * r2_x * r2**(-1.5)
    term2_y = (m1+m2)/r2 * (-time_transform(r1, r2, A1, A2)**2) * r2_y * r2**(-1.5) 
    term2_z = (m1+m2)/r2 * (-time_transform(r1, r2, A1, A2)**2) * r2_z * r2**(-1.5) 
    
    acceleration1j[2, 0] = - (term1_x + term2_x)
    acceleration1j[2, 1] = - (term1_y + term2_y)
    acceleration1j[2, 2] = - (term1_z + term2_z)
    
    # (2/4) Calculate k-th x, y, z accelerations from second part of the 
    #       Interaction Hamiltonian (in inertial coords)
    
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
                rk_minus_rj = np.sqrt(dx**2 + dy**2 + dz**2)
                rk_minus_rj3 = np.sqrt(dx**2 + dy**2 + dz**2)**3
                
                term1_x = time_transform(r1, r2, A1, A2) * mj * (rj_x-rk_x) / rk_minus_rj3
                term1_y = time_transform(r1, r2, A1, A2) * mj * (rj_y-rk_y) / rk_minus_rj3
                term1_z = time_transform(r1, r2, A1, A2) * mj * (rj_z-rk_z) / rk_minus_rj3
                term2_x = -mj/rk_minus_rj * (-time_transform(r1, r2, A1, A2)**2) * r2_x * r2**(-1.5)
                term2_y = -mj/rk_minus_rj * (-time_transform(r1, r2, A1, A2)**2) * r2_y * r2**(-1.5)
                term2_z = -mj/rk_minus_rj * (-time_transform(r1, r2, A1, A2)**2) * r2_z * r2**(-1.5)

                acceleration2[k, 0] += - (term1_x + term2_x)
                acceleration2[k, 1] += - (term1_y + term2_y)
                acceleration2[k, 2] += - (term1_z + term2_z)
            
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
                ri_minus_rk = np.sqrt(dx**2 + dy**2 + dz**2)
                ri_minus_rk3 = np.sqrt(dx**2 + dy**2 + dz**2)**3
                
                term1_x = time_transform(r1, r2, A1, A2) * mi * (rk_x-ri_x) / ri_minus_rk3
                term1_y = time_transform(r1, r2, A1, A2) * mi * (rk_y-ri_y) / ri_minus_rk3
                term1_z = time_transform(r1, r2, A1, A2) * mi * (rk_z-ri_z) / ri_minus_rk3
                term2_x = -mi/ri_minus_rk * (-time_transform(r1, r2, A1, A2)**2) * r2_x * r2**(-1.5)
                term2_y = -mi/ri_minus_rk * (-time_transform(r1, r2, A1, A2)**2) * r2_y * r2**(-1.5)
                term2_z = -mi/ri_minus_rk * (-time_transform(r1, r2, A1, A2)**2) * r2_z * r2**(-1.5)

                acceleration2[k, 0] -= - (term1_x + term2_x)
                acceleration2[k, 1] -= - (term1_y + term2_y)
                acceleration2[k, 2] -= - (term1_z + term2_z)
    
    # (3/4) Convert accelerations from second part of the Interaction Hamiltonian to Jacobi coords
    
    acceleration2j = vanillaWH_tools.inertial_to_jacobi_acc(sim, acceleration2)
    
    # (4/4) kick Jacobi velocities by h * (accel_part1 + accel_part2)
    
    for i in range(len(sim_jacobi)):
        
        sim_jacobi[i, 3] += h * (acceleration1j[i, 0] + acceleration2j[i, 0]) 
        sim_jacobi[i, 4] += h * (acceleration1j[i, 1] + acceleration2j[i, 1]) 
        sim_jacobi[i, 5] += h * (acceleration1j[i, 2] + acceleration2j[i, 2]) 
    
    return sim_jacobi