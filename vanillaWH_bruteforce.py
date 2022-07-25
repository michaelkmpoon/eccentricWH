"""vanillaWH_bruteforce.py: Wisdom-Holman integrator solving the EOM with an ODE solver"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint
import vanillaWH_tools

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
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(sim_jacobi[1,0]**2 + sim_jacobi[1,1]**2 + sim_jacobi[1,2]**2)
    r2 = np.sqrt(sim_jacobi[2,0]**2 + sim_jacobi[2,1]**2 + sim_jacobi[2,2]**2)
    
    mu1 = m1*m2 / (m1 + m2)
    mu2 = m3 * (m1 + m2) / (m1 + m2 + m3)
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    t = np.array([0, h])
    initial_vector1 = [r1_x, r1_y, r1_z, m2*v1_x, m2*v1_y, m2*v1_z]
    initial_vector2 = [r2_x, r2_y, r2_z, m3*v2_x, m3*v2_y, m3*v2_z]
    
    sol1 = odeint(drift_ODE1, initial_vector1, t, args=(r1, r2, m1, m2, m3, mu1, mu2))
    sol2 = odeint(drift_ODE2, initial_vector2, t, args=(r1, r2, m1, m2, m3, mu1, mu2))
    
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

def drift_ODE1(vector, t, r1, r2, m1, m2, m3, mu1, mu2):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
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
    
    rx_eqn = px/mu2
    ry_eqn = py/mu2
    rz_eqn = pz/mu2

    px_eqn = -m3*(m1+m2) * r2**(-3) * rx
    py_eqn = -m3*(m1+m2) * r2**(-3) * ry
    pz_eqn = -m3*(m1+m2) * r2**(-3) * rz
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt
