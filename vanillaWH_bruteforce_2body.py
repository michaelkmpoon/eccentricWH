"""vanillaWH_bruteforce.py: Wisdom-Holman integrator solving the EOM with an ODE solver"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import integrator_tools

def drift_solve_ivp(sim_jacobi, sim, h):
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

    r1_x, r1_y, r1_z, v1_x, v1_y, v1_z = sim_jacobi[1,:6]
    m1, m2 = sim[:2,6]

    total_mass = m1 + m2

    # Solve ODE for r' and p' in Jacobi coordinates

    t = np.array([0., h])
    initial_vector1 = [r1_x, r1_y, r1_z, v1_x, v1_y, v1_z]

    sol1 = solve_ivp(drift_ODE1_solve_ivp, t, initial_vector1, method='RK45', t_eval = t, rtol=1e-13, atol=1e-13, args=(total_mass,))

    # Update position and velocity with h*r' and h*p' in Jacobi coordinates

    r1_x = sol1.y[0,-1]
    r1_y = sol1.y[1,-1]
    r1_z = sol1.y[2,-1]
    v1_x = sol1.y[3,-1]
    v1_y = sol1.y[4,-1]
    v1_z = sol1.y[5,-1]

    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1_z, v1_x, v1_y, v1_z])

    return sim_jacobi

def drift_ODE1_solve_ivp(t, vector, total_mass):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """

    rx, ry, rz, vx, vy, vz = vector

    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(rx**2 + ry**2 + rz**2)

    rx_eqn = vx
    ry_eqn = vy
    rz_eqn = vz

    vx_eqn = -total_mass * r1**(-3.) * rx
    vy_eqn = -total_mass * r1**(-3.) * ry
    vz_eqn = -total_mass * r1**(-3.) * rz

    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn]

    return dvectordt

def drift_odeint(sim_jacobi, sim, h):
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
    m1, m2 = sim[:2,6]

    total_mass = m1 + m2
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    t = np.array([0., h])
    initial_vector1 = [r1_x, r1_y, r1_z, v1_x, v1_y, v1_z]
    
    sol1 = odeint(drift_ODE1_odeint, initial_vector1, t, rtol=1e-13, atol=1e-13, args=(total_mass,))

    # Update position and velocity with h*r' and h*p' in Jacobi coordinates
    
    r1_x = sol1[-1,0]
    r1_y = sol1[-1,1]
    r1_z = sol1[-1,2]
    v1_x = sol1[-1,3]
    v1_y = sol1[-1,4]
    v1_z = sol1[-1,5]
    
    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1_z, v1_x, v1_y, v1_z])
    
    return sim_jacobi

def drift_ODE1_odeint(vector, t, total_mass):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """

    rx, ry, rz, vx, vy, vz = vector
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(rx**2 + ry**2 + rz**2)
    
    rx_eqn = vx
    ry_eqn = vy
    rz_eqn = vz

    vx_eqn = -total_mass * r1**(-3.) * rx
    vy_eqn = -total_mass * r1**(-3.) * ry
    vz_eqn = -total_mass * r1**(-3.) * rz
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn]
    
    return dvectordt

def drift_euler(sim_jacobi, sim, h):
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
    m1, m2 = sim[:2,6]
    m1_jacobi, m2_jacobi = sim_jacobi[:2,6]
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(sim_jacobi[1,0]**2 + sim_jacobi[1,1]**2 + sim_jacobi[1,2]**2)
    
    #mu1 = m1*m2 / (m1 + m2)
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    #t = np.array([0., h])
    #initial_vector1 = [r1_x, r1_y, r1_z, m2_jacobi*v1_x, m2_jacobi*v1_y, m2_jacobi*v1_z]
    
    #sol1 = odeint(drift_ODE1, initial_vector1, t, args=(r1, m1, m2, mu1))

    # Update position and velocity with h*r' and h*p' in Jacobi coordinates
    
    r1_x += h*v1_x
    r1_y += h*v1_y
    r1_z += h*v1_z
    v1_x -= h*m1*r1**(-3.) * r1_x
    v1_y -= h*m1*r1**(-3.) * r1_y
    v1_z -= h*m1*r1**(-3.) * r1_z
    
    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1_z, v1_x, v1_y, v1_z])
    
    return sim_jacobi
