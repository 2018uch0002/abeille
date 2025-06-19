#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:20:54 2025

@author: singhp10
"""

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import h5py

# =============================================================================
# open the pre-calculated the file of the target tally and fine-mesh tally
# =============================================================================
file = h5py.File("c5g7_traget.h5", 'r')


x_bounds = np.array(file["x-bounds"])
fine_mesh_x_bounds = np.copy(x_bounds)

y_bounds = np.array(file["y-bounds"])
fine_mesh_y_bounds = np.copy(y_bounds)

NE_group = 7


# create the linear-array of tally-avg according to element-index
# get the number of element 
Nx = np.shape(x_bounds)[0] - 1
Ny = np.shape(y_bounds)[0] - 1

N_nodes = (Nx+1)*(Ny+1)

target_tally = np.array(file["c5g7-target-mesh"])

fine_mesh_tally_avg = np.array(file['fine-mesh'])

# =============================================================================
# get the 2D at the mid point of the fine mesh tally
# =============================================================================

fine_mesh_x = []
fine_mesh_y = []
dx_fine_mesh = fine_mesh_x_bounds[1] - fine_mesh_x_bounds[0]
for i in range(0, len(fine_mesh_x_bounds)-1):
    fine_mesh_x.append( 0.5 * (fine_mesh_x_bounds[i]+fine_mesh_x_bounds[i+1]) )

dy_fine_mesh = fine_mesh_y_bounds[1] - fine_mesh_y_bounds[0]
for i in range(0, len(fine_mesh_y_bounds)-1):
    fine_mesh_y.append( 0.5 * (fine_mesh_y_bounds[i]+fine_mesh_y_bounds[i+1]) )


# fine_mesh_avg_volume = dx_fine_mesh * dy_fine_mesh 

# fine_mesh_tally_avg /= fine_mesh_avg_volume


# =============================================================================
# plot the results
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


def evaluate_tally(g : int, xi: float, eta: float, ix : int, iy: int):
    value = 0.

    N1 = 0.25 * (1-xi) * (1-eta)
    value += target_tally[g, ix, iy, 0] * N1
    
    N2 = 0.25 * (1+xi) * (1-eta)
    value += target_tally[g, ix, iy, 1] * N2
    
    N3 = 0.25 * (1+xi) * (1+eta)
    value += target_tally[g, ix, iy, 2] * N3
    
    N4 = 0.25 * (1-xi) * (1+eta)
    value += target_tally[g, ix, iy, 3] * N4
    
    return value 


# get the points along the diagonal direction
diagonal_line = []
for ixy in range(0, len(fine_mesh_y)):
    point_y = fine_mesh_y[ixy]    
    point_x = fine_mesh_x[ixy]
    radial_value = np.sqrt(point_x** 2 + point_y**2)
    if (point_x < 0. and point_y < 0.):
        radial_value *= -1
    diagonal_line.append( radial_value )
        
        
# for g in range(0, NE_group):
#     # while reconstructing, flux needs to be normalised by the volume integral of the element    
#     linear_lagrange_reconstruct = np.zeros([len(fine_mesh_x), len(fine_mesh_y)])
    
#     for iy in range(0, len(fine_mesh_y)):
#         point_y = fine_mesh_y[iy]
                
#         for ix in range(0, len(fine_mesh_x)):
#             point_x = fine_mesh_x[ix]        
            
#             loc_ix, loc_iy = get_element_location_in_lagrange_linear(point_x, point_y)
            
#             xi, eta, area = get_scaled_xi_eta(point_x, loc_ix, point_y, loc_iy)
            
#             loc_ix, loc_iy = get_element_location_in_lagrange_linear(point_x, point_y)
            
#             linear_lagrange_reconstruct[ix, iy] = evaluate_tally(g, xi, eta, loc_ix, loc_iy) 
        
# # =============================================================================
# # plot the results
# # =============================================================================

#     plt.figure("lagrange-2D")
#     plt.pcolormesh(fine_mesh_x_bounds, fine_mesh_y_bounds, linear_lagrange_reconstruct)
    
#     plt.xlabel("x [cm]")
#     plt.ylabel("y [cm]")
#     plt.title("C5G7 group-{:} Lagrange Linear Quad4 Element".format(g))
    
#     plt.savefig("plots/flux_contour_group_{:}.png".format(g), dpi = 300, bbox_inches = "tight")
  

    
#     # plot along the diagonal line
#     flux_along_diag_lagrange = []
#     flux_along_fine_mesh = []
#     for ixy in range(0, len(fine_mesh_y)):
#         flux_along_diag_lagrange.append(linear_lagrange_reconstruct[ixy, ixy])    
#         flux_along_fine_mesh.append(fine_mesh_tally_avg[g, ixy, ixy])

#     plt.figure("line-plot-along-diagaonal-group-{:}".format(g))
#     plt.plot(diagonal_line, flux_along_diag_lagrange, label = "Lagrange Linear Quad4 Element")
#     plt.plot(diagonal_line, flux_along_fine_mesh, label = "standard-mesh-tally")
    
#     plt.xlabel("r [cm]")
#     plt.ylabel("Flux [Arb. Units]")
#     plt.title("C5G7: group-{:} flux along the diagonal".format(g))
#     plt.legend()
#     plt.show()
    # plt.savefig("plots/flux_along_diagonal_group_{:}.png".format(g), dpi = 300, bbox_inches = "tight")
      
  
# =============================================================================
# Plot the volume integral
# =============================================================================

# volume integral of N1, N2, N3, and N4 will be 1.
linear_lagrange_vol_integral = 0.25 * np.sum(target_tally, axis= 3) 
rel_error_linear_lagrange_vol_integral = 100. * (1 - np.divide(linear_lagrange_vol_integral, 
                                                               fine_mesh_tally_avg, 
                                                               out = np.zeros_like(linear_lagrange_vol_integral),
                                                               where = fine_mesh_tally_avg != 0.))

g = 0
index_truncate = 17*3*2 + 25*0
for g in range(0, NE_group):
    plt.figure("relative-error-2D")
    norm = TwoSlopeNorm(vcenter= 0., 
                        vmax = np.max(rel_error_linear_lagrange_vol_integral[g, : index_truncate, : index_truncate ]),
                        vmin = np.min(rel_error_linear_lagrange_vol_integral[g, : index_truncate, : index_truncate ])
                        )
    plt.pcolormesh(fine_mesh_x_bounds[: index_truncate+ 1], 
                   fine_mesh_y_bounds[: index_truncate + 1], 
                   rel_error_linear_lagrange_vol_integral[g, : index_truncate, : index_truncate],
                   cmap = "bwr",
                   norm = norm)
    plt.colorbar()
    
    plt.xlabel("x [cm]")
    plt.ylabel("y [cm]")
    plt.title("C5G7 group-{:} Relative Error (%)in the Lagrange Linear".format(g))
    
    # plt.savefig("plots/rel_error_lux_contour_group_{:}.png".format(g), dpi = 300, bbox_inches = "tight")
    plt.show()



  
  
