"""eccentricWH_bruteforce.py: Wisdom-Holman integrator (using ODE solver) for highly eccentric orbits following Mikkola 1997"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint
import integrator_tools

def drift(sim_jacobi, sim, A1, p0, h):
    """
    Advance (Keplerian-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        A1 (float): time transform parameters
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
        initial_vector = [rx, ry, rz, vx, vy, vz]
        sol = odeint(drift_ODE, initial_vector, t, rtol=1e-13, atol=1e-13, args=(total_mass, A1, p0))

        # Update position and velocity in Jacobi coordinates
        rx = sol[-1,0]
        ry = sol[-1,1]
        rz = sol[-1,2]
        vx = sol[-1,3]
        vy = sol[-1,4]
        vz = sol[-1,5]
        
        sim_jacobi[i,:6] = np.array([rx, ry, rz, vx, vy, vz])

    return sim_jacobi

def drift_ODE(vector, t, total_mass, A1, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, vx, vy, vz = vector
    
    # separation between the centre-of-mass of objects interior to the orbiting object, 
    # and the orbiting object
    separation = np.sqrt(rx**2 + ry**2 + rz**2)
    velocity_squared = vx**2 + vy**2 + vz**2
    
    time_transform_a1 = time_transform(separation, A1)
    rx_eqn = time_transform_a1*vx
    ry_eqn = time_transform_a1*vy
    rz_eqn = time_transform_a1*vz
    
    term1 = -time_transform_a1*total_mass * separation**(-3)
    term2 = - 0.5*velocity_squared + total_mass/separation - p0
    
    vx_eqn = term1*rx + time_transform_derivative(separation, A1, rx)*term2
    vy_eqn = term1*ry + time_transform_derivative(separation, A1, ry)*term2 
    vz_eqn = term1*rz + time_transform_derivative(separation, A1, rz)*term2 
          
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn]
    
    return dvectordt

# def drift_ODE(vector, t, total_mass, A1, p0):
#     """
#     EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
#     according to eqn. (53) of Mikkola 1997.
#     This is for the first non-central object.
#     """
    
#     r, v = vector[:3], vector[3:]
    
#     # separation between the centre-of-mass of objects interior to the orbiting object, 
#     # and the orbiting object
#     separation = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
#     velocity_squared = v[0]**2 + v[1]**2 + v[2]**2
    
#     time_transform_a1 = time_transform(separation, A1)
#     r_eqn = [time_transform_a1 * v_component for v_component in v]
    
#     term1 = -time_transform_a1*total_mass * separation**(-3)
#     term2 = - 0.5*velocity_squared + total_mass/separation - p0
    
#     v_eqn = [term1*r_component + time_transform_derivative(separation, A1, r_component)*term2 for r_component in r]
    
#     dvectordt = r_eqn + v_eqn
    
#     return dvectordt

def time_transform(separation, A1):
    return separation/A1
    #return 1.

def time_transform_derivative(separation, A1, ri):    
    #return (1/A1) * (ri/separation)
    return 0.
