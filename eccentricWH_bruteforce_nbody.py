"""eccentricWH_bruteforce.py: Wisdom-Holman integrator (using ODE solver) for highly eccentric orbits following Mikkola 1997"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint
import integrator_tools

def drift(sim_jacobi, sim, A0, A1, A2, p0, h):
    """
    Advance (Keplerian-like part) system by timestep h,
    according to eqn. (53) of Mikkola 1997.
    
    Parameters:
        sim_jacobi (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        A0, A1, A2 (float): time transform parameters
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
        
        if i == 1:
            r2 = np.sqrt(sim_jacobi[2,0]**2+sim_jacobi[2,1]**2+sim_jacobi[2,2]**2)
            args = (total_mass, A0, A1, A2, p0, i, r2)
        elif i == 2:
            r1 = np.sqrt(sim_jacobi[1,0]**2+sim_jacobi[1,1]**2+sim_jacobi[1,2]**2)
            args = (total_mass, A0, A1, A2, p0, i, r1)
        else:
            r1 = np.sqrt(sim_jacobi[1,0]**2+sim_jacobi[1,1]**2+sim_jacobi[1,2]**2)
            r2 = np.sqrt(sim_jacobi[2,0]**2+sim_jacobi[2,1]**2+sim_jacobi[2,2]**2)
            args = (total_mass, A0, A1, A2, p0, i, np.array([r1, r2]))
            
        sol = odeint(drift_ODE, initial_vector, t, rtol=1e-13, atol=1e-13, args=args)
            
        # Update position and velocity in Jacobi coordinates
        rx = sol[-1,0]
        ry = sol[-1,1]
        rz = sol[-1,2]
        vx = sol[-1,3]
        vy = sol[-1,4]
        vz = sol[-1,5]
        
        sim_jacobi[i,:6] = np.array([rx, ry, rz, vx, vy, vz])

    return sim_jacobi

def drift_ODE(vector, t, total_mass, A0, A1, A2, p0, object_num, constant_r):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    #This is for the first non-central object.
    """
    
    rx, ry, rz, vx, vy, vz = vector
    
    # separation between the centre-of-mass of objects interior to the orbiting object, 
    # and the orbiting object
    separation = np.sqrt(rx**2 + ry**2 + rz**2)
    velocity_squared = vx**2 + vy**2 + vz**2
    
    if object_num == 1:
        r1 = separation
        r2 = constant_r
    elif object_num == 2:
        r1 = constant_r
        r2 = separation
    else:
        r1 = constant_r[0]
        r2 = constant_r[1]

    time_transformed = time_transform(r1, r2, A0, A1, A2)
    rx_eqn = time_transformed*vx
    ry_eqn = time_transformed*vy
    rz_eqn = time_transformed*vz
    
    term1 = - time_transformed * total_mass * separation**(-3)
    term2 = - 0.5*velocity_squared + total_mass/separation - p0[object_num-1]
    #print(object_num, - 0.5*velocity_squared, total_mass/separation, -p0[object_num-1])
    
    vx_eqn = term1*rx + time_transform_derivative(r1, r2, A0, A1, A2, object_num, rx)*term2
    vy_eqn = term1*ry + time_transform_derivative(r1, r2, A0, A1, A2, object_num, ry)*term2 
    vz_eqn = term1*rz + time_transform_derivative(r1, r2, A0, A1, A2, object_num, rz)*term2 
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, vx_eqn, vy_eqn, vz_eqn]
    
    return dvectordt

def time_transform(r1, r2, A0, A1, A2):
    return 1/(A0 + A1/r1 + A2/r2)

def time_transform_derivative(r1, r2, A0, A1, A2, object_num, ri): 
    if object_num == 1:
        return (A0 + A1/r1 + A2/r2)**(-2.) * (A1*ri) / r1**3
    elif object_num == 2:
        return (A0 + A1/r1 + A2/r2)**(-2.) * (A2*ri) / r2**3
    else:
        return 0.
    #return 0.
