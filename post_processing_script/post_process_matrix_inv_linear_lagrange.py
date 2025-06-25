"""
This post-process script to solve the mass-matrix constructed from the 
elemental matrix to reconstruct the flux
"""

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import h5py

# =============================================================================
# open the miel and coeff of lagrange linear at the nodes points
# =============================================================================
file = h5py.File("miel.h5", 'r')

lagrange_tally = file["results"]["lagrange-linear"]

position_filter_id = lagrange_tally.attrs["position-filter"]
position_filter = file["tally-filters"]["position-filters"][str(position_filter_id)]
x_bounds = np.array(position_filter["x-bounds"])
y_bounds = np.array(position_filter["y-bounds"])

NE_group = 1
is_energy_filter = False
if ("energy-filter" in lagrange_tally.attrs.keys()):
    is_energy_filter = True
    energy_filter_id = lagrange_tally.attrs["energy-filter"]
    energy_filter = file["tally-filters"]["energy-filters"][str(energy_filter_id)]
    energy_bounds = np.array(energy_filter["energy-bounds"])
    NE_group = np.shape(energy_bounds)[0]-1

tally_avg = np.array(lagrange_tally["avg"])

# create the linear-array of tally-avg according to element-index
# get the number of element 
Nx = np.shape(x_bounds)[0] - 1
Ny = np.shape(y_bounds)[0] - 1

tally_avg_shape = np.array([NE_group, Nx, Ny, 4])

N_nodes = (Nx+1)*(Ny+1)
source_tally = np.zeros([NE_group, N_nodes])
          
for g in range(0, NE_group):
    for ix in range(0, Nx+1):
        for iy in range(0, Ny+1):
            i_node = ix + iy * (Nx+1)
            if (is_energy_filter == True):
                source_tally[g, i_node] = tally_avg[g, ix, iy]
            else: 
                source_tally[g, i_node] = tally_avg[ix, iy]



"""# the below method shall be helpful, if in the monte carlo simulation
# tally is being evaluated for the nodes in a element, but tally is 
# not beind done on shared nodes unlike above method. Here, each 
# element's nodes will kept isolated while tallying, and then summed
# in the post process.  
for g in range(0, NE_group):
    for ix in range(0, Nx):
        for iy in range(0, Ny):
            if (is_energy_filter == True):
                i_node = ix + iy * (Nx+1)
                source_tally[g, i_node] += tally_avg[g, ix, iy, 0] # N1
                
                i_node = (ix+1) + iy * (Nx+1)
                source_tally[g, i_node] += tally_avg[g, ix, iy, 1] # N2
                
                i_node = (ix+1) + (iy+1) * (Nx+1)
                source_tally[g, i_node] += tally_avg[g, ix, iy, 2] # N3
    
                i_node = ix + (iy+1) * (Nx+1)
                source_tally[g, i_node] += tally_avg[g, ix, iy, 3] # N4
            else:
                i_node = ix + iy * (Nx+1)
                source_tally[g, i_node] += tally_avg[ix, iy, 0] # N1
                
                i_node = (ix+1) + iy * (Nx+1)
                source_tally[g, i_node] += tally_avg[ix, iy, 1] # N2
                
                i_node = (ix+1) + (iy+1) * (Nx+1)
                source_tally[g, i_node] += tally_avg[ix, iy, 2] # N3
    
                i_node = ix + (iy+1) * (Nx+1)
                source_tally[g, i_node] += tally_avg[ix, iy, 3] # N4
"""

# =============================================================================
# construction of global mass-matrix
# =============================================================================

# elemental mass matrix
M_element = np.array([[4, 2, 1, 2],
                      [2, 4, 2, 1],
                      [1, 2, 4, 2],
                      [2, 1, 2, 4]])

"""
6 | 7 | 8
3 | 4 | 5 
0 | 1 | 2

ix >>  and  iy ^

linear_index = 

"""
global_mass_matrix = np.zeros([N_nodes, N_nodes])

# go over each element
for iy in range(0, Ny):
    for ix in range(0, Nx):
        # area of the element
        Ae = (x_bounds[ix+1] - x_bounds[ix])*(y_bounds[iy+1] - y_bounds[iy])
        
        # node at xmin-ymin of the element
        N1 = ix + iy * (Nx+1)
        
        # node at xmax-ymin of the element
        N2 = (ix+1) + iy * (Nx+1)
        
        # node at xmax-ymax of the element
        N3 = (ix+1) + (iy+1) * (Nx+1)
        
        # node at xmin-ymax of the element
        N4 = ix + (iy+1) * (Nx+1)
        
        node_index = [N1, N2, N3, N4]
        ele_row = 0
        for N_node_x in [N1, N2, N3, N4]:
            ele_column = 0
            for N_node_y in [N1, N2, N3, N4]:
                
                global_mass_matrix[N_node_x, N_node_y] += M_element[ele_row, ele_column] *(Ae/36)
                ele_column += 1
            ele_row += 1        

inv_global_mass_matrix = np.linalg.inv(global_mass_matrix)

# =============================================================================
# Get the solution of X (which are at nodes) by solving the MX=b or X = inv(M)B 
# =============================================================================

target_tally = np.zeros(tally_avg_shape)
for g in range(0, NE_group):
    target_tally_vector = np.dot(inv_global_mass_matrix, source_tally[g, :])
    
    # the target tally is arragened element wise so, reconstruct it to 2D array
    # also it needs to be normalised by the volume integral
    for ix in range(0, Nx):
        for iy in range(0, Ny):
            # N1
            i_node = ix + iy * (Nx+1)
            target_tally[g, ix, iy, 0] = target_tally_vector[i_node]
            
            # N2
            i_node = (ix+1) + iy * (Nx+1)
            target_tally[g, ix, iy, 1] = target_tally_vector[i_node]
            
            # N3
            i_node = (ix+1) + (iy+1) * (Nx+1)
            target_tally[g, ix, iy, 2] = target_tally_vector[i_node]
    
            # N4
            i_node = ix + (iy+1) * (Nx+1)
            target_tally[g, ix, iy, 3] = target_tally_vector[i_node]

# store the target_tally
file_target_tally = h5py.File("miel_lagrange_tally.h5", 'w')
file_target_tally.create_dataset("elemental-coefficients", data=target_tally)
file_target_tally.create_dataset("x-bounds", data = x_bounds)
file_target_tally.create_dataset("y-bounds", data = y_bounds)
if (is_energy_filter):
    file_target_tally.create_dataset("energy-bounds", data = energy_bounds)
file_target_tally.close()

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
