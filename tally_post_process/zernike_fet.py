#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 09:03:00 2024

@author: singhp10
"""

import h5py
import numpy as np
import filter as flt
from zernike_polynomial import ZernikePolynomials as ZrPoly
from scipy.special import eval_legendre

# =============================================================================
# Zernike-FET Reconstruction
# =============================================================================
class ZernikeFET:
    def __init__(self, filename: str ="", tally_name: str ="", file : h5py.File = ''):
        
        if ( file != "" and isinstance(file, h5py.File) ):
            self.file = file
        elif ( filename != "" and isinstance(filename, str) ):
            self.__filename = filename
            self.__file = h5py.File(filename, "r")
        
        else:
            raise TypeError("incorrect file or filename is given.")
        
        if ( tally_name == "" or not isinstance(tally_name, str) ):
            raise TypeError("invalid tally-name is given.")
        
        self.__tally_name = tally_name
        if ( self.__tally_name not in self.__file['results'].keys() ):
            raise TypeError("given tally-name \""  + tally_name + "\" is not found.")
        
        self.__tally = self.__file['results'][self.__tally_name]
        self.__tally_avg = self.__file['results'][self.__tally_name]['avg']
        self.__tally_std = self.__file['results'][self.__tally_name]['std']
        
        # now get the tally-type
        tally_type_list = np.array(self.__file['results'][self.__tally_name].attrs['type'])
        self.type = tally_type_list[tally_type_list != 0].tobytes().decode("ascii")
        
        # get the poisiton and the energy-filter if id exist
        if ( 'energy-filter' in self.__file['results'][self.__tally_name].attrs.keys() ):
            energy_filter_id = self.__file['results'][self.__tally_name].attrs['energy-filter']
            energy_bounds = list( self.__file['tally-filters']['energy-filters'][str(energy_filter_id)]['energy-bounds'] )
            self.__energy_filter = flt.EnergyFilter(energy_bounds, energy_filter_id)
        else: 
            self.__energy_filter = None
            
        if ( 'position-filter' in self.__file['results'][self.__tally_name].attrs.keys() ):
            position_filter_id = self.__file['results'][self.__tally_name].attrs['position-filter']
            pos_ftr = self.__file['tally-filters']['position-filters'][str(position_filter_id)].attrs
            name_array = np.array( self.__file['tally-filters']['position-filters'][str(position_filter_id)].attrs['type'] )
            self.position_filter_type = name_array[ name_array != 0 ].tobytes().decode('ascii')
        else:
            raise TypeError("position-filter is not found for the tally-type " + self.type)
            
        if ( self.position_filter_type == 'cylinder-filter' ):
            name_array = np.array(pos_ftr['axis'])
            axis = name_array[ name_array != 0 ].tobytes().decode('ascii')
            origin = pos_ftr['origin']
            origin_point = flt.Position(origin[0], origin[1], origin[2])
            pitch = pos_ftr['pitch']
            radius = pos_ftr['radius']
            shape = pos_ftr['shape']
            self.__position_filter = flt.CylinderFilter(origin_point, radius, pitch[0], pitch[1], pitch[2], shape[0], shape[1], shape[2], axial_axis= axis , filter_id = position_filter_id)
        else:
            raise TypeError("incorrect position-filter is given. it must be a cylinder-filter")
        
        self.zernike_order = int(self.__tally.attrs['zernike-order'])
        if ( 'legendre-order' in self.__tally.attrs.keys() ):
            self.legendre_order = int(self.__tally.attrs['legendre-order'])
        else:
            self.legendre_order = None
        
        self.zr_polynomial = ZrPoly(self.zernike_order)
        
        
    def evaluate_zernike_fet(self, energy_and_position : list, is_divide_volum = True) -> list:
        
        value = []
        i = -1
        for r_E in energy_and_position:
            i += 1
            E = r_E[0]
            r = flt.Position(r_E[1], r_E[2], r_E[3])
            tally_avg = self.__tally_avg
            # get the enrgy-bin if exists.
            if ( self.__energy_filter != None ):
                E_index = self.__energy_filter.get_index(E)
                if ( E_index[0] == True ):
                    tally_avg = tally_avg[E_index[1]]
                    
                else:
                    value.append(0.)
                    continue
            
            # get the position-index
            pos_index = self.__position_filter.get_index(r)
            if ( pos_index[0] == False ):
                value.append(0.)
                continue

            for i in pos_index[1]:
                tally_avg = tally_avg[i]
            
            (scaled_r, theta) = self.__position_filter.get_scaled_radius_and_center(pos_index[1], r)
            print(i, pos_index[1], scaled_r, r)
            zr_poly_values = self.zr_polynomial.eval_zernike(scaled_r, theta)
            tally_zr = 0.
            for order in range(0, self.zernike_order+1):
                tally_zr += tally_avg[0][order] * zr_poly_values[order] / self.zr_polynomial.orthonormalisation_constant(order)
            
            tally_legendre = 1.
            inv_phi0 = 1.
            if ( self.legendre_order != None ):
                zmin = self.__position_filter.z_min(pos_index[1])
                zmax = self.__position_filter.z_max(pos_index[1])
                scaled_z = 2. * ( r.z() - zmin ) / ( zmax - zmin )
                tally_legendre = 0.
                for order in range(0, self.legendre_order+1):
                    tally_legendre += (2*order+1) * tally_avg[1][order] * eval_legendre(order, scaled_z)
                
                inv_phi0 = 1. / tally_avg[0][0]
                            
            inv_dV = 1.
            if ( is_divide_volum == True ):
                inv_dV = 1. / self.__position_filter.dV()

            value.append( tally_zr * tally_legendre * inv_dV * inv_phi0)

        return value    
            
            
            
        