"""
This script can be only run once the lobal elemental mass matrix 
has been inverted and had the process the actual coefficent of
linear shape function on Quad4 element. 
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# =============================================================================
# open the pre-calculated the file of the target tally and fine-mesh tally
# =============================================================================
file_lagrange = h5py.File("miel_lagrange_tally.h5", 'r')

x_bounds = np.array(file_lagrange["x-bounds"])
y_bounds = np.array(file_lagrange["y-bounds"])

NE_group = 1
is_energy_filter = False
if ("energy-bounds" in file_lagrange.keys()):
    is_energy_filter = True
    energy_bounds = file_lagrange["energy-bounds"]


# create the linear-array of tally-avg according to element-index
# get the number of element 
Nx = np.shape(x_bounds)[0] - 1
Ny = np.shape(y_bounds)[0] - 1

N_nodes = (Nx+1)*(Ny+1)

target_tally = np.array(file_lagrange["elemental-coefficients"])

file_lagrange.close()


# =============================================================================
# Functions to evaluate the corrdinate and Lagrange Tally
# =============================================================================

# function to get the scaled x and y or xi and eta
def get_scaled_xi_eta(x: float, loc_ix : int, 
                      y: float, loc_iy : int ):
    x_ele_bounds = [x_bounds[loc_ix], x_bounds[loc_ix+1]]
    y_ele_bounds = [y_bounds[loc_iy], y_bounds[loc_iy+1]]

    xi = 2. * (x - x_ele_bounds[0]) / (x_ele_bounds[1] - x_ele_bounds[0]) - 1
    eta = 2. * (y - y_ele_bounds[0]) / (y_ele_bounds[1] - y_ele_bounds[0]) - 1.
    
    area = (x_ele_bounds[1] - x_ele_bounds[0]) * (y_ele_bounds[1] - y_ele_bounds[0])

    return xi, eta, area

dx_lagrange_linear = x_bounds[1] - x_bounds[0]
dy_lagrange_linear = y_bounds[1] - y_bounds[0]

def get_element_location_in_lagrange_linear(x: float, y:float):
    ix = int( (x-x_bounds[0]) / dx_lagrange_linear )
    iy = int( (y-y_bounds[0]) / dy_lagrange_linear )    
    
    return (ix, iy)

def evaluate_tally(g : int, point_x: float, point_y: float):
    # get the index 
    loc_ix, loc_iy = get_element_location_in_lagrange_linear(point_x, point_y)
    
    # get the scaled position
    xi, eta, area = get_scaled_xi_eta(point_x, loc_ix, point_y, loc_iy)
    
    value = 0.

    N1 = 0.25 * (1-xi) * (1-eta)
    value += target_tally[g, loc_ix, loc_iy, 0] * N1
    
    N2 = 0.25 * (1+xi) * (1-eta)
    value += target_tally[g, loc_ix, loc_iy, 1] * N2
    
    N3 = 0.25 * (1+xi) * (1+eta)
    value += target_tally[g, loc_ix, loc_iy, 2] * N3
    
    N4 = 0.25 * (1-xi) * (1+eta)
    value += target_tally[g, loc_ix, loc_iy, 3] * N4
    
    return value 


# =============================================================================
# get the 2D at the mid point of the fine mesh tally and comapres 
# =============================================================================

file = h5py.File("miel.h5", 'r')
fine_mesh_tally = file["results"]["fine-mesh"]
fine_mesh_tally_avg = np.array(fine_mesh_tally["avg"])
fine_mesh_tally_pos_filter_id = fine_mesh_tally.attrs["position-filter"]
fine_mesh_position_filter = file["tally-filters"]["position-filters"][str(fine_mesh_tally_pos_filter_id)]
fine_mesh_x_bounds = fine_mesh_position_filter["x-bounds"]
fine_mesh_y_bounds = fine_mesh_position_filter["y-bounds"]
fine_mesh_z_bounds = fine_mesh_position_filter["z-bounds"]

fine_mesh_x = []
fine_mesh_y = []
dx_fine_mesh = fine_mesh_x_bounds[1] - fine_mesh_x_bounds[0]
for i in range(0, len(fine_mesh_x_bounds)-1):
    fine_mesh_x.append( 0.5 * (fine_mesh_x_bounds[i]+fine_mesh_x_bounds[i+1]) )

dy_fine_mesh = fine_mesh_y_bounds[1] - fine_mesh_y_bounds[0]
for i in range(0, len(fine_mesh_y_bounds)-1):
    fine_mesh_y.append( 0.5 * (fine_mesh_y_bounds[i]+fine_mesh_y_bounds[i+1]) )

dz_fine_mesh = fine_mesh_z_bounds[1] - fine_mesh_z_bounds[0]

fine_mesh_avg_volume = dx_fine_mesh * dy_fine_mesh #* dz_fine_mesh
fine_mesh_tally_avg /= fine_mesh_avg_volume

for g in range(0, NE_group):
    linear_lagrange_reconstruct = np.zeros([len(fine_mesh_x), len(fine_mesh_y)])
    
    for iy in range(0, len(fine_mesh_y)):
        point_y = fine_mesh_y[iy]
                
        for ix in range(0, len(fine_mesh_x)):
            point_x = fine_mesh_x[ix]        
            linear_lagrange_reconstruct[ix, iy] = evaluate_tally(g, point_x, point_y) 
        

# =============================================================================
# Relative Difference in Volume integral of the quantity
# =============================================================================

# # volume integral of N1, N2, N3, and N4 will be 1.
# linear_lagrange_vol_integral = 0.25 * np.sum(target_tally, axis= 3) 
# rel_error_linear_lagrange_vol_integral = 100. * (1 - np.divide(linear_lagrange_vol_integral, 
#                                                             fine_mesh_tally_avg, 
#                                                             out = np.zeros_like(linear_lagrange_vol_integral),
#                                                             where = fine_mesh_tally_avg != 0.))
