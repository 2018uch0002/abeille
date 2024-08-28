#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 01:24:10 2024

@author: parthsingh
"""
import numpy as np
import h5py
import filter as flt
from scipy.special import eval_legendre

class TallyPostProcess:
    def __init__(self, filename : str = "" , tally_name: str = "",  file: h5py.File = ""):
        
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
        
        elif ( self.type != 'general-tally' ):
            raise TypeError("position-filter is not found for the tally-type " + self.type)
        
        else: 
            self.__energy_filter = None
        
        # construct the poistion-filter accordig to the type
        if ( self.position_filter_type == 'regular-cartesian-mesh' ):
            low = pos_ftr['low']
            low_point = flt.Position(low[0], low[1], low[2])
            high = pos_ftr['high']
            high_point = flt.Position(high[0], high[1], high[2])
            shape = pos_ftr['shape']
            Nx = shape[0]
            Ny = shape[1]
            Nz = shape[2]
            self.__position_filter = flt.RegularReatilearFilter(low_point, high_point, Nx, Ny, Nz, position_filter_id)
        
        elif ( self.position_filter_type == 'cylinder-filter' ):
            name_array = np.array(pos_ftr['axis'])
            axis = name_array[ name_array != 0 ].tobytes().decode('ascii')
            origin = pos_ftr['origin']
            origin_point = flt.Position(origin[0], origin[1], origin[2])
            pitch = pos_ftr['pitch']
            radius = pos_ftr['radius']
            shape = pos_ftr['shape']
            self.__position_filter = flt.CylinderFilter(origin_point, radius, pitch[0], pitch[1], pitch[2], shape[0], shape[1], shape[2], axial_axis= axis , filter_id = position_filter_id)

    def evaluate( self, energy_and_position : list, is_divide_volum = True) -> list:
        """
        \'evaluate\' funtion can evlauate the tally if a list of energy and position combined by a tuple of shape 4, where 
        shape[0] : energy
        shape[1] : x
        shape[2] : y
        shape[3] : z
        
        is provided. The function will retrun a list 
        """
        if self.type == 'general-tally' :
            return self.general_tally_evaluate( energy_and_position, is_divide_volum )
        
        if self.type == "legendre-fet" :
            return self.legendre_fet_evaluate(energy_and_position, is_divide_volum )
        
    
    def general_tally_evaluate(self, energy_and_position : list , is_divide_volum = True) -> list:
        """
        \'general_tally_evaluate\' funtion can evlauate the tally if a list of energy and position combined by a tuple of shape 4, where 
        shape[0] : energy
        shape[1] : x
        shape[2] : y
        shape[3] : z
        is provided. The function will retrun a list 
        """
        value = []

        for r_E in energy_and_position:

            E = r_E[0]
            r = flt.Position(r_E[1], r_E[2], r_E[3])
            index = []
            tally_avg = self.__tally_avg
            # get the energy_index if exists            
            if ( self.__energy_filter != None ):
                E_index = self.__energy_filter.get_index(E)
                if (E_index[0] == True):
                    tally_avg = tally_avg[E_index[1]]
                else:
                    value.append(0.)
                    continue
            
            # get the poisiton index if exist
            if (self.__position_filter != None):
                pos_index = self.__position_filter.get_index(r)
                if ( pos_index[0] == True ):
                    for i in pos_index[1]:
                        tally_avg = tally_avg[i]
                    
                    if (is_divide_volum == True):
                        inv_dV = self.__position_filter.inv_dV()
                        tally_avg *= inv_dV
                
                    value.append(tally_avg)
                    
                else:
                    value.append(0)
                    continue
                
        return value
                    
    def legendre_fet_evaluate(self, energy_and_position : list , is_divide_volum = True) -> list:
        """
        \'legendre_fet_evaluate\' funtion can evlauate the tally if a list of energy and position combined by a tuple of shape 4, where 
        shape[0] : energy
        shape[1] : x
        shape[2] : y
        shape[3] : z
        is provided. The function will retrun a list 
        """
        value = []
        axes = self.__tally.attrs['axes']
        legendre_orders = list(self.__tally.attrs['order'])
        
        for r_E in energy_and_position:
            E = r_E[0]
            r = flt.Position(r_E[1], r_E[2], r_E[3])
            index = []
            tally_avg = self.__tally_avg
            # get the energy_index if exists            
            if ( self.__energy_filter != None ):
                E_index = self.__energy_filter.get_index(E)
                if (E_index[0] == True):
                    tally_avg = tally_avg[E_index[1]]
                else:
                    value.append(0.)
                    continue
            

            pos_index = self.__position_filter.get_index(r)
            if ( pos_index[0] == True ):
                for i in pos_index[1]:
                    tally_avg = tally_avg[i]
            else:
                value.append(0)
            
                continue
            #print(pos_index, "\t", end="")
            # itrate over the axis and order
            tally_value = 1.
            it_coeff = 0;
            phi0 = tally_avg[0]
            for it_axis in range(0, len(axes)):
                
                if ( axes[it_axis] == 'x' or axes[it_axis] == "X" ):
                    xmin = self.__position_filter.x_min(pos_index[1])
                    inv_dx = self.__position_filter.inv_dx()
                    scaled_x = (r.x() - xmin) * inv_dx * 2. - 1.
                    fet_x_value = 0
                    # print(xmin)
                    for order in range(0, int(legendre_orders[it_axis])+1):
                        fet_x_value += ( 2*order + 1 ) * tally_avg[it_coeff] * eval_legendre(order, scaled_x)
                        it_coeff += 1
                    tally_value *= fet_x_value

                
                if ( axes[it_axis] == 'y' or axes[it_axis] == "Y" ):
                    ymin = self.__position_filter.y_min(pos_index[1])
                    inv_dy = self.__position_filter.inv_dy()
                    scaled_y = (r.y() - ymin) * inv_dy * 2. - 1.
                    fet_y_value = 0
                    for order in range(0, int(legendre_orders[it_axis])+1):
                        fet_y_value += ( 2*order + 1 ) * tally_avg[it_coeff] * eval_legendre(order, scaled_y)
                        it_coeff += 1
                    tally_value *= fet_y_value

                    
                if ( axes[it_axis] == 'z' or axes[it_axis] == "Z" or axes[it_axis] == 's'):
                    zmin = self.__position_filter.z_min(pos_index[1])
                    inv_dz = self.__position_filter.inv_dz()
                    scaled_z = (r.z() - zmin) * inv_dz * 2. - 1.
                    fet_z_value = 0
                    for order in range(0, int(legendre_orders[it_axis])+1):
                        fet_z_value += ( 2*order + 1 ) * tally_avg[it_coeff] * eval_legendre(order, scaled_z)
                        it_coeff += 1
                    tally_value *= fet_z_value

                    
                if ( it_axis > 0 ):
                    tally_value /= phi0
            inv_dV = self.__position_filter.inv_dV()
            tally_value *= inv_dV
            value.append(tally_value)
            
        return value
        
        
    def close(self):
        self.__file.close()
    
        
    def evaluate_std( self, energy_and_position : list, is_divide_volum = True) -> list:
        """
        \'evaluate\' funtion can evlauate the tally if a list of energy and position combined by a tuple of shape 4, where 
        shape[0] : energy
        shape[1] : x
        shape[2] : y
        shape[3] : z
        
        is provided. The function will retrun a list 
        """
        if self.type == 'general-tally' :
            return self.general_tally_evaluate_std( energy_and_position, is_divide_volum )
        
    
    def general_tally_evaluate_std(self, energy_and_position : list , is_divide_volum = True) -> list:
        """
        \'general_tally_evaluate\' funtion can evlauate the tally if a list of energy and position combined by a tuple of shape 4, where 
        shape[0] : energy
        shape[1] : x
        shape[2] : y
        shape[3] : z
        is provided. The function will retrun a list 
        """
        value = []
        for r_E in energy_and_position:
            E = r_E[0]
            r = flt.Position(r_E[1], r_E[2], r_E[3])
            tally_std = self.__tally_std
            # get the energy_index if exists            
            if ( self.__energy_filter != None ):
                E_index = self.__energy_filter.get_index(E)
                if (E_index[0] == True):
                    tally_std = tally_std[E_index[1]]
                else:
                    value.append(0.)
                    continue
            
            # get the poisiton index if exist
            if (self.__position_filter != None):
                pos_index = self.__position_filter.get_index(r)
                if ( pos_index[0] == True ):
                    for i in pos_index[1]:
                        tally_std = tally_std[i]
                    
                    if (is_divide_volum == True):
                        inv_dV = self.__position_filter.inv_dV()
                        tally_std *= inv_dV
                
                    value.append(tally_std)
                    
                else:
                    value.append(0)
                    continue
                
        return value    