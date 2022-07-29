"""eccentricWH_bruteforce.py: Wisdom-Holman integrator (using ODE solver) for highly eccentric orbits following Mikkola 1997"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import solve_ivp
import integrator_tools

def drift(sim_jacobi, sim, A1, physical_time, p0, h):
    """
    Advance (Keplerian-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        A1 (float): time transform parameters
        physical_time (float): time coordinate q0 as described in Mikkola 1997
        p0 (float): negative of Hamiltonian at time 0
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
        initial_vector = [rx, ry, rz, vx, vy, vz, physical_time, p0]
        sol = solve_ivp(drift_ODE, t, initial_vector, method='RK45', t_eval = t, rtol=1e-13, atol=1e-13, args=(total_mass, A1, physical_time, p0))

        # Update position and velocity in Jacobi coordinates
        rx = sol.y[0,-1]
        ry = sol.y[1,-1]
        rz = sol.y[2,-1]
        vx = sol.y[3,-1]
        vy = sol.y[4,-1]
        vz = sol.y[5,-1]
        physical_time = sol.y[6,-1]
        
        sim_jacobi[i,:6] = np.array([rx, ry, rz, vx, vy, vz])
    
    return sim_jacobi, physical_time

def drift_ODE(t, vector, total_mass, A1, physical_time, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, vx, vy, vz, physical_time, p0 = vector
    
    # separation between the centre-of-mass of objects interior to the orbiting object, 
    # and the orbiting object
    separation = np.sqrt(rx**2 + ry**2 + rz**2)
    velocity_squared = vx**2 + vy**2 + vz**2
    
    rx_eqn = time_transform(separation, A1)*vx
    ry_eqn = time_transform(separation, A1)*vy
    rz_eqn = time_transform(separation, A1)*vz
    
    term1 = -time_transform(separation, A1)*total_mass * separation**(-3)
    term2 = - 0.5*velocity_squared + total_mass/separation - p0
    vx_eqn = term1*rx - time_transform_derivative(separation, A1)*term2
    vy_eqn = term1*ry - time_transform_derivative(separation, A1)*term2 
    vz_eqn = term1*rz - time_transform_derivative(separation, A1)*term2
    q0_eqn = time_transform(separation, A1)
    p0_eqn = 0.
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn, q0_eqn, p0_eqn]
    
    return dvectordt

def time_transform(separation, A1):
    return separation/A1
    #return 1.

def time_transform_derivative(separation, A1):
    return 1/A1
    #return 0.
