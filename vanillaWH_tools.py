"""vanillaWH_tools.py: move to centre-of-mass, coordinate transforms, calculate energy"""

__author__ = 'Michael Poon'
__email__ = 'michaelkm.poon@mail.utoronto.ca'

import numpy as np

def move_to_com(sim):
    """
    Convert from heliocentric to inertial (centre-of-mass) positions and velocities
    
    Parameters:
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in heliocentric coordinates
        
    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coodinates
    """    
    
    
    m_total = 0.
    COM_x, COM_y, COM_z, COM_vx, COM_vy, COM_vz = np.zeros(6)
    
    for i in range(len(sim)):
        
        xi, yi, zi, vxi, vyi, vzi, mi = sim[i]
        
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
    
    r0_x, r0_y, r0_z, v0_x, v0_y, v0_z, m0 = sim[0]
    
    COM_x = m0 * r0_x
    COM_y = m0 * r0_y
    COM_z = m0 * r0_z
    COM_vx = m0 * v0_x
    COM_vy = m0 * v0_y
    COM_vz = m0 * v0_z
    
    M_iminus1 = 0.

    for i in range(1,len(sim)):
        
        ri_x, ri_y, ri_z, vi_x, vi_y, vi_z, mi = sim[i]      

        # x, y, z, vx, vy, vz in jacobi coordinates for the ith object 
        
        m_iminus1 = sim[i-1, 6]
        M_iminus1 += m_iminus1
        
        rij_x = ri_x - COM_x/M_iminus1
        rij_y = ri_y - COM_y/M_iminus1
        rij_z = ri_z - COM_z/M_iminus1
        vij_x = vi_x - COM_vx/M_iminus1
        vij_y = vi_y - COM_vy/M_iminus1
        vij_z = vi_z - COM_vz/M_iminus1
        
        simj[i,:6] = rij_x, rij_y, rij_z, vij_x, vij_y, vij_z
        
        COM_x = COM_x * (1 + mi/M_iminus1) + mi*rij_x
        COM_y = COM_y * (1 + mi/M_iminus1) + mi*rij_y
        COM_z = COM_z * (1 + mi/M_iminus1) + mi*rij_z
        COM_vx = COM_vx * (1 + mi/M_iminus1) + mi*vij_x
        COM_vy = COM_vy * (1 + mi/M_iminus1) + mi*vij_y
        COM_vz = COM_vz * (1 + mi/M_iminus1) + mi*vij_z
        
        # calculate Jacobi mass
        mij = (mi*M_iminus1) / (mi+M_iminus1)
        simj[i, 6] = mij
        
    # Different convention for 0-th pos/vel
    m_total = M_iminus1 + mi
    r0j_x = COM_x / m_total
    r0j_y = COM_y / m_total
    r0j_z = COM_z / m_total
    v0j_x = COM_vx / m_total
    v0j_y = COM_vy / m_total
    v0j_z = COM_vz / m_total

    
    # Jacobi mass of central object is total mass
    m0j = 0.
    for i in range(len(sim)):
        mi = sim[i, 6]
        m0j += mi
       
    simj[0] = r0j_x, r0j_y, r0j_z, v0j_x, v0j_y, v0j_z, m0j
    
    return simj

def inertial_to_jacobi_acc(sim, accel_list):
    """
    Convert from inertial (centre-of-mass) accelerations to Jacobi accelerations
    
    This follows transformations.c from Rebound, WHFast paper, and Mikkola & Innanen 1999
    
    Parameters:
        sim (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
        acceleration_list (array): list of accelerations from 2nd part of Interaction Hamiltonian
        
    Returns:
        (array): list of Jacobi accelerations
    """
    
    accelerationj = np.zeros((len(accel_list), 3))
    
    m0 = sim[0, 6]
    a0_x, a0_y, a0_z = accel_list[0]
    
    COM_ax = m0 * a0_x
    COM_ay = m0 * a0_y
    COM_az = m0 * a0_z

    M_iminus1 = 0.

    for i in range(1,len(sim)):
        
        mi = sim[i, 6]
        
        m_iminus1 = sim[i-1, 6]
        M_iminus1 += m_iminus1
        
        ai_x, ai_y, ai_z = accel_list[i]

        # ax, ay, az in jacobi coordinates for the ith object 
        aij_x = ai_x - COM_ax/M_iminus1
        aij_y = ai_y - COM_ay/M_iminus1
        aij_z = ai_z - COM_az/M_iminus1
        
        accelerationj[i] = np.array([aij_x, aij_y, aij_z])
        
        COM_ax = COM_ax * (1 + mi/M_iminus1) + mi*aij_x
        COM_ay = COM_ay * (1 + mi/M_iminus1) + mi*aij_y
        COM_az = COM_az * (1 + mi/M_iminus1) + mi*aij_z
    
    # Different convention for 0-th coordinate
    m_total = M_iminus1 + mi
    a0j_x = COM_ax / m_total
    a0j_y = COM_ay / m_total
    a0j_z = COM_az / m_total
    
    accelerationj[0] = np.array([a0j_x, a0j_y, a0j_z])
    
    return accelerationj

def jacobi_to_inertial(simj, masses):
    """
    Convert from Jacobi positions and velocities to inertial pos/vel
    
    Parameters:
        simj (array): [x, y, z, vx, vy, vz, m] stacked for each particle in Jacobi coordinates
        
    Returns:
        (array): [x, y, z, vx, vy, vz, m] stacked for each particle in inertial coordinates
    """
    
    sim = np.zeros(np.shape(simj))
    
    m0 = masses[0]
    r0j_x, r0j_y, r0j_z, v0j_x, v0j_y, v0j_z = simj[0,:6]

    m_total = 0.
    for i in range(len(sim)):
        mi = masses[i]
        m_total += mi

    COM_x = r0j_x * m_total
    COM_y = r0j_y * m_total
    COM_z = r0j_z * m_total
    COM_vx = v0j_x * m_total
    COM_vy = v0j_y * m_total
    COM_vz = v0j_z * m_total
    
    Mi = 0.
    
    for i in range(len(simj)-1, 0, -1):
        
        mi = masses[i]
        Mi = np.sum(masses[:i+1])
        M_iminus1 = Mi - mi
        
        rij_x, rij_y, rij_z, vij_x, vij_y, vij_z = simj[i,:6]
        
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
        
        sim[i,:6] = np.array([ri_x, ri_y, ri_z, vi_x, vi_y, vi_z])
        
        COM_x *= M_iminus1
        COM_y *= M_iminus1
        COM_z *= M_iminus1
        COM_vx *= M_iminus1
        COM_vy *= M_iminus1
        COM_vz *= M_iminus1
    
    # Set the inertial pos/vel for the central object
    r0_x = COM_x / m0
    r0_y = COM_y / m0
    r0_z = COM_z / m0
    v0_x = COM_vx / m0
    v0_y = COM_vy / m0
    v0_z = COM_vz / m0

    sim[0,:6] = np.array([r0_x, r0_y, r0_z, v0_x, v0_y, v0_z])
    sim[:,6] = masses
    
    return sim

def energy_fn(sim):
    """ Calculate energy of system as shown in p. 77 of Tremaine's book
    
    Parameters:
        sim (rebound.simulation.Simulation): 
            Rebound simulation object

    Returns:
        (float): Total energy of the system 
    """
   
    e_kin = 0.
    e_pot = 0.
    
    for i in range(len(sim)):
        
        xi, yi, zi, vxi, vyi, vzi, mi = sim[i]
        
        e_kin += 0.5 * mi * (vxi**2 + vyi**2 + vzi**2)
        
        for j in range(len(sim)):
            
            xj, yj, zj, mj = sim[j,0], sim[j,1], sim[j,2], sim[j,6]
        
            if i > j:
                dx = xi - xj
                dy = yi - yj
                dz = zi - zj
                e_pot += - mi * mj / np.sqrt(dx*dx + dy*dy + dz*dz)
    
    return e_kin + e_pot
