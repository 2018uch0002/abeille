#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 12:37:15 2025

@author: singhp10
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
file = h5py.File("miel.h5", 'r')

fine_mesh_tally = file["results"]["mesh-tally"]
fine_mesh_tally_avg = np.array(fine_mesh_tally["avg"])
fine_mesh_tally_pos_filter_id = fine_mesh_tally.attrs["position-filter"]
fine_mesh_position_filter = file["tally-filters"]["position-filters"][str(fine_mesh_tally_pos_filter_id)]
fine_mesh_x_bounds = fine_mesh_position_filter["x-bounds"]
fine_mesh_y_bounds = fine_mesh_position_filter["y-bounds"]
fine_mesh_z_bounds = fine_mesh_position_filter["z-bounds"]

g = 6
plt.pcolormesh(fine_mesh_x_bounds, fine_mesh_y_bounds, fine_mesh_tally_avg[g, :, :])
plt.xlabel("x [cm]")
plt.ylabel("y [cm]")

plt.title("group-{:}".format(g))
plt.savefig("standard-mesh-tally-group-{:}.png".format(g), dpi = 300, bbox_inches = "tight")