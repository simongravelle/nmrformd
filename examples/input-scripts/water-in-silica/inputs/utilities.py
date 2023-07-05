import copy
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R

def generate_random_location(Lx, Ly, Lz, recenter=True):
    """
    generate a random location within a given box
    """  
    if recenter == True:
        x = np.random.rand()*Lx-Lx/2
        y = np.random.rand()*Ly-Ly/2
        z = np.random.rand()*Lz-Lz/2
    else:
        x = np.random.rand()*Lx
        y = np.random.rand()*Ly
        z = np.random.rand()*Lz   
    return x, y, z

def generate_random_orientation(XYZ_initial):
    """
    generate 3D aleatory rotation for molecule coordinate
    """
    rotation_degrees = np.random.rand()*90
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([np.random.rand(), 
                              np.random.rand(), 
                              np.random.rand()])
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    XYZ_rotated = rotation.apply(XYZ_initial)
    return XYZ_rotated


def search_closest_neighbor(XYZ_neighbor, XYZ_molecule, Lx, Ly, Lz):
    """
        Search neighbor in a box and return the closest distance with a molecule
        
        Periodic boundary conditions are automaticaly accounted
    """
    box = np.array([Lx, Ly, Lz])
    min_distance = np.max(box)
    for XYZ_atom in XYZ_molecule:
        dxdydz = np.remainder(XYZ_neighbor - XYZ_atom + box/2., box) - box/2.
        min_distance = np.min([min_distance,np.min(norm(dxdydz,axis=1))])
    return min_distance
