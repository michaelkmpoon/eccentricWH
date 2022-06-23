import rebound 
# import matplotlib.pyplot as plt
import numpy as np
import time
import math
from scipy import optimize

def drift(simj, sim, object_num, h):
    """
    Advance (Keplerian part) system by timestep h
    
    Parameters:
        simj (array): [x, y, vx, vy, m] stacked for each particle in Jacobi coordinates
        sim (array): [x, y, vx, vy, m] stacked for each particle in inertial coordinates
        object_num (int): index of orbiting particle
        h (float): timestep

    Returns:
        (array): [x, y, vx, vy, m] stacked for each updated particle in Jacobi coodinates  
    """
    
    r = simj[object_num, :3]
    v = simj[object_num, 3:6]
    M = 0.
    for i in range(object_num+1):
        M += sim[i, 6]
    
    r0 = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    r0i = 1./r0  #r0 inverse
    v2 = v[0]**2 + v[1]**2 + v[2]**2
    beta = 2*M*r0i - v2
    eta0 = r[0]*v[0] + r[1]*v[1] + r[2]*v[2]
    zeta0 = M - beta*r0

    # Solve Kepler's equation for X
    dtr0i = h * r0i # first order guess
    X_guess = dtr0i * (1 - dtr0i*eta0*0.5*r0i) # second order guess (according to Rebound)
    
    X = optimize.newton(func=Keplers_eqn, x0=X_guess, args=(r0, beta, eta0, zeta0, h))
    
    G1 = calculate_G(1, beta, X) 
    G2 = calculate_G(2, beta, X)
    G3 = calculate_G(3, beta, X)
    
    ri = 1./(r0 + eta0*G1 + zeta0*G2)
    f = 1 - M*G2*r0i
    g = h - M*G3
    fd = -M*G1*r0i*ri   # time derivative of f
    gd = 1 - M*G2*ri    # time derivative of g
    
    # solve for new position and velocity using the f and g functions
    rx = f*r[0] + g*v[0]
    ry = f*r[1] + g*v[1]
    rz = f*r[2] + g*v[2]
    vx = fd*r[0] + gd*v[0]
    vy = fd*r[1] + gd*v[1]
    vz = fd*r[2] + gd*v[2]

    simj[object_num, 0] = rx
    simj[object_num, 1] = ry
    simj[object_num, 2] = rz
    simj[object_num, 3] = vx
    simj[object_num, 4] = vy
    simj[object_num, 5] = vz

    return simj

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

def move_to_com(sim):
    """
    Convert from heliocentric to inertial (centre-of-mass) positions and velocities
    
    Parameters:
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in heliocentric coordinates
        
    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coodinates
    """
    
    m_total = 0.
    COM_x, COM_y, COM_z = 0., 0., 0.
    COM_vx, COM_vy, COM_vz = 0., 0., 0.
    
    for i in range(len(sim)):
        
        mi = sim[i, 6]
        xi, yi, zi = sim[i, 0], sim[i, 1], sim[i, 2]
        vxi, vyi, vzi = sim[i, 3], sim[i, 4], sim[i, 5]
        
        m_total += mi
        COM_x += mi*xi
        COM_y += mi*yi
        COM_z += mi*zi
        COM_vx += mi*vxi
        COM_vy += mi*vyi
        COM_vz += mi*vzi
        
    COM_x /= m_total
    COM_y /= m_total
    COM_z /= m_total
    COM_vx /= m_total
    COM_vy /= m_total
    COM_vz /= m_total
    
    for i in range(len(sim)):
        sim[i, 0] -= COM_x
        sim[i, 1] -= COM_y
        sim[i, 2] -= COM_z
        sim[i, 3] -= COM_vx
        sim[i, 4] -= COM_vy
        sim[i, 5] -= COM_vz
        
    return sim
    
def inertial_to_jacobi(sim):
    """
    Convert from inertial (centre-of-mass) positions and velocities to Jacobi pos/vel
    
    Jacobi Coordinates: coordinates are measured relative to the COM of all inner bodies
    
    This follows transformations.c from Rebound, WHFast paper, and Mikkola & Innanen 1999
    
    Parameters:
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        
    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coodinates
    """
    
    simj = np.zeros(np.shape(sim))
    
    m0 = sim[0, 6]
    r0 = sim[0, :3]
    v0 = sim[0, 3:6]
    
    COM_x = m0 * r0[0]
    COM_y = m0 * r0[1]
    COM_z = m0 * r0[2]
    COM_vx = m0 * v0[0]
    COM_vy = m0 * v0[1]
    COM_vz = m0 * v0[2]
    
    M_iminus1 = 0.

    for i in range(1,len(sim)):
        
        mi = sim[i, 6]
        
        m_iminus1 = sim[i-1, 6]
        M_iminus1 += m_iminus1
        
        ri_x = sim[i, 0]
        ri_y = sim[i, 1]
        ri_z = sim[i, 2]
        vi_x = sim[i, 3]
        vi_y = sim[i, 4]
        vi_z = sim[i, 5]
        
        # x, y, vx, vy in jacobi coordinates for the ith object 
        rij_x = ri_x - COM_x/M_iminus1
        rij_y = ri_y - COM_y/M_iminus1
        rij_z = ri_z - COM_z/M_iminus1
        vij_x = vi_x - COM_vx/M_iminus1
        vij_y = vi_y - COM_vy/M_iminus1
        vij_z = vi_z - COM_vz/M_iminus1
        
        simj[i, 0] = rij_x
        simj[i, 1] = rij_y
        simj[i, 2] = rij_z
        simj[i, 3] = vij_x
        simj[i, 4] = vij_y
        simj[i, 5] = vij_z
        
        COM_x = COM_x * (1 + mi/M_iminus1) + mi*rij_x
        COM_y = COM_y * (1 + mi/M_iminus1) + mi*rij_y
        COM_z = COM_z * (1 + mi/M_iminus1) + mi*rij_z
        COM_vx = COM_vx * (1 + mi/M_iminus1) + mi*vij_x
        COM_vy = COM_vy * (1 + mi/M_iminus1) + mi*vij_y
        COM_vz = COM_vz * (1 + mi/M_iminus1) + mi*vij_z
        
        # Calculate Jacobi mass
        simj[i, 6] = (mi*M_iminus1) / (mi+M_iminus1)
    
    # Different convention for 0-th coordinate
    m_total = M_iminus1 + mi
    r0j_x = COM_x / m_total
    r0j_y = COM_y / m_total
    r0j_z = COM_z / m_total
    v0j_x = COM_vx / m_total
    v0j_y = COM_vy / m_total
    v0j_z = COM_vz / m_total
    
    simj[0, 0] = r0j_x
    simj[0, 1] = r0j_y
    simj[0, 2] = r0j_z
    simj[0, 3] = v0j_x
    simj[0, 4] = v0j_y
    simj[0, 5] = v0j_z
    
    # Jacobi mass of central object is total mass
    m0j = 0.
    for i in range(len(sim)):
        mi = sim[i, 6]
        m0j += mi
        
    simj[0, 6] = m0j
    
    return simj

def jacobi_to_inertial(simj, sim):
    """
    Convert from Jacobi positions and velocities to inertial pos/vel
    
    Parameters:
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        
    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coodinates
    """
    
    m0 = sim[0, 6]
    r0j = simj[0, :3]
    v0j = simj[0, 3:6]
    
    m_total = 0.
    for i in range(len(sim)):
        mi = sim[i, 6]
        m_total += mi

    COM_x = r0j[0] * m_total
    COM_y = r0j[1] * m_total
    COM_z = r0j[2] * m_total
    COM_vx = v0j[0] * m_total
    COM_vy = v0j[1] * m_total
    COM_vz = v0j[2] * m_total
    
    Mi = 0.
    
    # Iterate over objects from outer planet to inner planet with a for loop
    for i in range(len(simj)-1, 0, -1):
        
        mi = sim[i, 6]
        Mi = np.sum(sim[:i+1, 6])
        M_iminus1 = Mi - mi
        
        rij_x = simj[i, 0]
        rij_y = simj[i, 1]
        rij_z = simj[i, 2]
        vij_x = simj[i, 3]
        vij_y = simj[i, 4]
        vij_z = simj[i, 5]
        
        COM_x = (COM_x - mi * rij_x) / Mi
        COM_y = (COM_y - mi * rij_y) / Mi
        COM_z = (COM_z - mi * rij_z) / Mi
        COM_vx = (COM_vx - mi * vij_x) / Mi
        COM_vy = (COM_vy - mi * vij_y) / Mi
        COM_vz = (COM_vz - mi * vij_z) / Mi
        
        ri_x = rij_x + COM_x
        ri_y = rij_y + COM_y
        ri_z = rij_z + COM_z
        vi_x = vij_x + COM_vx
        vi_y = vij_y + COM_vy
        vi_z = vij_z + COM_vz
        
        sim[i, 0] = ri_x
        sim[i, 1] = ri_y
        sim[i, 2] = ri_z
        sim[i, 3] = vi_x
        sim[i, 4] = vi_y
        sim[i, 5] = vi_z
        
        COM_x *= M_iminus1
        COM_y *= M_iminus1
        COM_z *= M_iminus1
        COM_vx *= M_iminus1
        COM_vy *= M_iminus1
        COM_vz *= M_iminus1
    
    # Set the inertial coordinate for the central object
    r0_x = COM_x / m0
    r0_y = COM_y / m0
    r0_z = COM_z / m0
    v0_x = COM_vx / m0
    v0_y = COM_vy / m0
    v0_z = COM_vz / m0

    sim[0, 0] = r0_x
    sim[0, 1] = r0_y
    sim[0, 2] = r0_z
    sim[0, 3] = v0_x
    sim[0, 4] = v0_y
    sim[0, 5] = v0_z
    
    return sim

def energy_fn(sim):
    """ Calculate energy per unit mass as shown in p. 77 of Tremaine's book
    
    Parameters:
        sim (rebound.simulation.Simulation): 
            Rebound simulation object

    Returns:
        (float): Energy per unit mass
    """
   
    e_kin = 0.
    e_pot = 0.
    
    for i in range(len(sim)):
        
        mi = sim[i,6]
        xi, yi, zi = sim[i,0], sim[i,1], sim[i,2]
        vxi, vyi, vzi = sim[i,3], sim[i,4], sim[i,5]
        e_kin += 0.5 * m * (vxi**2 + vyi**2 + vzi**2)
        
        for j in range(len(sim)):
            
            xj, yj, zj, mj = sim[j,0], sim[j,1], sim[j,2], sim[j,6]
        
            if i > j:
                dx = xi - xj
                dy = yi - yj
                dz = zi - zj
                e_pot += - sim.G * m / np.sqrt(dx*dx + dy*dy + dz*dz)
    
    return e_kin + e_pot
