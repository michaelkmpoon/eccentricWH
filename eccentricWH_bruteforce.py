"""eccentricWH.py: Wisdom-Holman integrator for highly eccentric orbits following Mikkola 1997"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np
import math
from scipy import optimize
from scipy.integrate import odeint

def drift(sim_jacobi, sim, p0, h):
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
    
    r1_x, r1_y, r1z, v1x, v1y, v1z = sim_jacobi[1,:6]
    r2_x, r2_y, r2z, v2x, v2y, v2z = sim_jacobi[2,:6]
    m1, m2, m3 = sim[:3,6]
    
    # Norm of Jacobi coordinates. will be used in time transform
    r1 = np.sqrt(sim_jacobi[1,0]**2 + sim_jacobi[1,1]**2 + sim_jacobi[1,2]**2)
    r2 = np.sqrt(sim_jacobi[2,0]**2 + sim_jacobi[2,1]**2 + sim_jacobi[2,2]**2)
    
    mu1 = m1*m2 / (m1 + m2)
    mu2 = m3 * (m1 + m2) / (m1 + m2 + m3)
    
    # Solve ODE for r' and p' in Jacobi coordinates
    
    sol1 = odeint(drift_ODE1, initial_vector, t, args=(r1, r2, m1, m2, mu1, p0))
    sol2 = odeint(drift_ODE2, initial_vector, t, args=(r1, r2, m1, m2, m3, mu2, p0))
    
    # Update position and velocity with h*r' and h*p' in Jacobi coordinates
    
    r1_x += h*sol1[:,0]
    r1_y += h*sol1[:,1]
    r1_z += h*sol1[:,2]
    v1_x += h*sol1[:,3]
    v1_y += h*sol1[:,4]
    v1_z += h*sol1[:,5]
    
    r2_x += h*sol2[:,0]
    r2_y += h*sol2[:,1]
    r2_z += h*sol2[:,2]
    v2_x += h*sol2[:,3]
    v2_y += h*sol2[:,4]
    v2_z += h*sol2[:,5]
    
    sim_jacobi[1,:6] = np.array([r1_x, r1_y, r1z, v1x, v1y, v1z])
    sim_jacobi[2,:6] = np.array([r2_x, r2_y, r2z, v2x, v2y, v2z])
    
    return sim_jacobi

def drift_ODE1(vector, t, r1, r2, m1, m2, mu1, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the first non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    rx_eqn = time_transform(r1, r2)*px/mu1
    ry_eqn = time_transform(r1, r2)*py/mu1
    rz_eqn = time_transform(r1, r2)*pz/mu1
    
    p_squared = (px**2 + py**2 + pz**2)
    Kepler_term1 = p_squared/(2*mu1) - m1*m2/r1
    Kepler_term2 = p_squared/(2*mu2) - m3*(m1+m2) / r2
    Kepler_terms = Kepler_term1 + Kepler_term2 + p0
    
    dKepler_terms_wrt_r1x = (m1*m2/r1) * r1**(-3) * rx
    dKepler_terms_wrt_r1y = (m1*m2/r1) * r1**(-3) * ry
    dKepler_terms_wrt_r1z = (m1*m2/r1) * r1**(-3) * rz
    
    term2 = Kepler_terms * time_transform_derivative1(r1, r2)
    
    px_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r1x + term2)
    py_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r1y + term2)
    pz_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r1z + term2)
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt

def drift_ODE2(vector, t, r1, r2, m1, m2, m3, mu2, p0):
    """
    EOM from Hamilton's equations for the first part of the new Hamiltonian (Gamma_0),
    according to eqn. (53) of Mikkola 1997.
    This is for the second non-central object.
    """
    
    rx, ry, rz, px, py, pz = vector
    
    rx_eqn = time_transform(r1, r2)*px/mu2
    ry_eqn = time_transform(r1, r2)*py/mu2
    rz_eqn = time_transform(r1, r2)*pz/mu2
    
    p_squared = (px**2 + py**2 + pz**2)
    Kepler_term1 = p_squared/(2*mu1) - m1*m2/r1
    Kepler_term2 = p_squared/(2*mu2) - m3*(m1+m2) / r2
    Kepler_terms = Kepler_term1 + Kepler_term2 + p0
    
    dKepler_terms_wrt_r2x = (m3*(m1+m2) / r2) * r1**(-3) * rx
    dKepler_terms_wrt_r2y = (m3*(m1+m2) / r2) * r1**(-3) * ry
    dKepler_terms_wrt_r2z = (m3*(m1+m2) / r2) * r1**(-3) * rz
    
    term2 = Kepler_terms * time_transform_derivative2(r1, r2)
    
    px_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r2x + term2)
    py_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r2y + term2)
    pz_eqn = - (time_transform(r1, r2) * dKepler_terms_wrt_r2z + term2)
    
    dvectordt = [rx_eqn, ry_eqn, rz_eqn, px_eqn, py_eqn, pz_eqn]
    
    return dvectordt

def time_transform(r1, r2):
    return 1. / (1./r1 + 15./r2) 

def time_transform_derivative1(r1, r2):
    """
    Derivative of time transform w.r.t. r1
    """
    return (1./r1 + 15./r2)**(-2) * r1**(-2) 

def time_transform_derivative2(r1, r2):
    """
    Derivative of time transform w.r.t. r2
    """
    return 15. * (1./r1 + 15./r2)**(-2) * r2**(-2) 

def kick(sim_jacobi, sim, p0, h):
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

    # (1/4) Calculate x, y, z accelerations from first part of the 
    #       Interaction-like Hamiltonian (in Jacobi coords)
    
    acceleration1j = np.zeros((len(sim_jacobi), 3))
    
    
    
    # (2/4) Calculate k-th x, y, z accelerations from second part of the 
    #       Interaction Hamiltonian (in inertial coords)
    
    # (3/4) Convert accelerations from second part of the Interaction Hamiltonian to Jacobi coords
    
    # (4/4) kick Jacobi velocities by h * (accel_part1 + accel_part2)
    
    return sim_jacobi