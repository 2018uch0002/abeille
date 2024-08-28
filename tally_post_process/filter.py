#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 08:26:19 2024

@author: parthsingh
"""
from enum import Enum
from typing import Union
import numpy as np
import math 


# =============================================================================
# Position (x, y, z)
# =============================================================================
class Position:
    def __init__(self, x: float, y: float, z: float):
        self.__x = x
        self.__y = y
        self.__z = z
    
    def __str__(self):
        return "{:.3f}".format(self.__x) + ", \t" + "{:.3f}".format(self.__y) + ", \t" + "{:.3f}".format(self.__z)
        
    def x(self) -> float :
        return self.__x
    def y(self) -> float :
        return self.__y
    def z(self) -> float :
        return self.__z


# =============================================================================
# Energy Filter 
# =============================================================================
class EnergyFilter:
    def __init__(self, energy_bounds: Union[list, np.ndarray] , id_num: int = 0 ):
        self.__energy_bounds = energy_bounds.copy()
        # check the energy_bounds is a correct instance or not
        if not ( isinstance(energy_bounds, list) or isinstance(energy_bounds, np.ndarray) ):
            raise TypeError(" for energy-filter, the energy-bounds must be array.")
            
        self.__id = id_num
        
        # check the energy-filter making condition 
        if ( len(self.__energy_bounds) < 1 ):
            raise TypeError("the length of the energy-bounds are less than 2 for energy-filer with id " + str(self.__id) + ".")
        if ( self.__energy_bounds != sorted(self.__energy_bounds) ):
            raise TypeError("the energy_bounds are not sorted for the energy-filter  with id " + str(self.__id) + ".")
        if ( self. __energy_bounds[0] < 0.):
            raise TypeError("the first element of energy-filter with id " + str(self.__id) + " is negative.")
        
    def type_str(self):
        return  "energy-filter"
        
    def id(self):
        return self.__id
    def size(self):
        return len(self.__energy_bounds)
    def get_energy_bounds(self):
        return self.__energy_bounds.copy()
    
    # get the tuple with first element True or False if element found
    # and second for index in the energy-bounds
    def get_index(self, E):
        for i in range(0, len(self.__energy_bounds)-1):
            if ( self.__energy_bounds[i] <= E and E < self.__energy_bounds[i+1] ):
                return (True, i)
            
        return (False, -1)


# =============================================================================
# Regular Reactilinear Cartesian Filter
# =============================================================================
class RegularReatilearFilter:
    def __init__(self, low_point: Position, high_point: Position, Nx: int, Ny: int, Nz: int, filter_id: int = 0):
        self.__id = filter_id
        self.__low = low_point
        # check weather the low_point is a correct instance or not
        if not isinstance(low_point, Position):
            raise TypeError("\"low\" should be the object of \"Position\" for the Regular Reactlinear Filter with id" + str(filter_id) + ".")
            
        self.__high = high_point
        # check weather the high_point is a correct instance or not        
        if not isinstance(high_point, Position):
            raise TypeError("\"high\" should be the object of \"Position\" for the Regular Reactlinear Filter with id" + str(filter_id) + ".")
            
        # check weather Nx, Ny or Nz are greater than zero or not
        if ( Nx < 0 or Ny < 0 or Nz < 0):
            raise TypeError("Number of bins are negative for the Regular Reactlinear Filter with id" + str(filter_id) + ".")
        self.__Nx = Nx
        self.__Ny = Ny
        self.__Nz = Nz
        self.__dx = (self.__high.x() - self.__low.x()) / self.__Nx
        self.__inv_dx = 1.0 / self.__dx
        self.__dy = (self.__high.y() - self.__low.y()) / self.__Ny
        self.__inv_dy = 1.0 / self.__dy
        self.__dz = (self.__high.z() - self.__low.z()) / self.__Nz
        self.__inv_dz = 1.0 / self.__dz
        
        self.__x_index = 0
        self.__y_index = 1
        self.__z_index = 2
        if ( self.__Nx == 1 ):
            self.__x_index = 0
            self.__y_index -= 1
            self.__z_index -= 1
        
        if ( self.__Ny == 1 ):
            self.__y_index = 0
            self.__z_index -= 1
        
        if ( self.__Nz == 1 ):
            self.__z_index = 0
            
        
    # get the index for the bin
    def get_index(self, point: Position):
        ix = math.floor( (point.x() - self.__low.x()) * self.__inv_dx )
        iy = math.floor( (point.y() - self.__low.y()) * self.__inv_dy )
        iz = math.floor( (point.z() - self.__low.z()) * self.__inv_dz )
        
        if ( ix >= 0 and ix < self.__Nx 
            and iy >= 0 and iy < self.__Ny
            and iz >= 0 and iz < self.__Nz ):
            
            return (True, self.__reduce_dimension(ix, iy, iz))
        
        return (False, [-1, -1, -1])
    
    def __reduce_dimension(self, ix: int, iy : int, iz :int):
        if ( self.__Nx == 1 and self.__Ny and self.__Nz ):
            return [ix]
        indices = []
        
        if ( self.__Nx > 1 ):
            indices.append(ix)
        
        if ( self.__Ny > 1 ):
            indices.append(iy)
        
        if ( self.__Nz > 1 ):
            indices.append(iz)
        
        return indices
    
    def x_min(self, indices: list):
        if (self.__Nx == 1):
            return self.__low.x()
        
        return self.__low.x() + indices[self.__x_index] * self.__dx

    def x_max(self, indices: list):
        if (self.__Nx == 1):
            return self.__high.x()
        
        return self.__low.x() + indices[self.__x_index] * self.__dx + self.__dx

    def y_min(self, indices: list):
        if (self.__Ny == 1):
            return self.__low.y()
        
        return self.__low.y() + indices[self.__y_index] * self.__dy

    def y_max(self, indices: list):
        if (self.__Ny == 1):
            return self.__high.y()
        
            return self.__low.y() + indices[self.__x_index] * self.__dy + self.__dy

    def z_min(self, indices: list):
        if (self.__Nz == 1):
            return self.__low.z()
        
        return self.__low.z() + indices[self.__z_index] * self.__dz

    def z_max(self, indices: list):
        if (self.__Nz == 1):
            return self.__high.z()
        
        return self.__low.z() + indices[self.__z_index] * self.__dz + self.__dz
    
    def type_str(self):
        return "regular-reactilinear-cartesian"
    
    def id(self):
        return self.__id
    
    def dx(self):
        return self.__dx
    
    def inv_dx(self):
        return self.__inv_dx
    
    def dy(self):
        return self.__dy
    
    def inv_dy(self):
        return self.__inv_dy
    
    def dz(self):
        return self.__dz
    
    def inv_dz(self):
        return self.__inv_dz
    
    def dV(self):
        return self.__dx * self.__dy * self.__dz

    def inv_dV(self):
        return self.__inv_dx * self.__inv_dy * self.__inv_dz

# =============================================================================
# Cylinder Filter
# =============================================================================
class CylinderFilter:
    # a private class for the orientation for the axial direction
    class __Orientation:
        X = False
        Y = False
        Z = False
        pass
    
    def __map_coordinate(self, point: Position):
        if ( self.__Orientation.Z == True ):
            return point
        
        elif ( self.__Orientation.X == True ):
            p = Position(point.z(), point.y(), point.x() )
            return p
        
        elif ( self.__Orientation().Y == True ):
            p = Position(point.x(), point.z(), point. y() )
            return p
    
    def __map_index(self, indices: list):
        if ( self.__Orientation.Z == True ):
            return indices
        
        elif ( self.__Orientation.X == True ):
            return [ indices[2], indices[1], indices[0] ]
        
        elif ( self.__Orientation.Y == True ):
            return [ indices[0], indices[2], indices[1] ]
        
    def __init__(self, origin_point: Position, radius: float,
                 pitch_x: float, pitch_y: float, pitch_z: float, 
                 Nx: int, Ny: int, Nz: int, axial_axis: str = 'Z', filter_id : int = 0):
        self.__id = filter_id
                
        # get the radius
        if ( radius < 0.):
            raise ValueError("radius must be positive for the cylinder filter with id " + str(filter_id) + ".")
        self.__radius = radius
            
        # check the pitches
        if ( pitch_x < 0 or pitch_y < 0 or pitch_z < 0 ):
            raise ValueError("the pitches must be >= 0. for cylinder filter with id " + str(filter_id) + ".")

        # get the axial-axis
        if ( axial_axis == "X" or axial_axis == "x" ):
            self.__Orientation.X = True
            self.__pitch_x = pitch_z
            self.__pitch_y = pitch_y
            self.__pitch_z = pitch_x
        elif ( axial_axis == "Y" or axial_axis == "y"):
            self.__Orientation.Y = True
            self.__pitch_x = pitch_x
            self.__pitch_y = pitch_z
            self.__pitch_z = pitch_y
        elif ( axial_axis == "Z" or axial_axis == "z" ):
            self.__Orientation.Z = True
            self.__pitch_x = pitch_x
            self.__pitch_y = pitch_y
            self.__pitch_z = pitch_z
        else:
            raise TypeError("invalid \"axis\"")
        
        self.__axis = axial_axis
        
        # check for the correct input for origin and map it
        if not isinstance(origin_point, Position):
            raise TypeError("origin shpuld be the object of \"Position\" for Cylinder-Filter with id " + str(filter_id) + ".")
        
        self.__origin = self.__map_coordinate(origin_point)
        
        # get the number of bins
        self.__Real_Nx = Nx
        self.__Real_Ny = Ny
        self.__Real_Nz = Nz
        
        # get the N, Ny and Nz oriented according to the cylinder orientation
        maped_shape = self.__map_index([Nx, Ny, Nz])
        self.__Nx = maped_shape[0]
        self.__Ny = maped_shape[1]
        self.__Nz = maped_shape[2]
        
        # check the bins condition and get the pitches 
        self.__infinte_length = False
        if ( Nx == 0 and self.__Orientation.X == True ):
            self.__infinte_length = True
        elif ( Nx <= 0 ):
            raise TypeError("the number of bin for x direction is not valid for cylinder filter with id " + str(filter_id) + ".")
        
        if ( Ny == 0 and self.__Orientation.Y == True ):
            self.__infinte_length = True
        elif ( Ny <= 0 ):
            raise TypeError("the number of bin for y direction is not valid for cylinder filter with id " + str(filter_id) + ".")
        
        if ( Nz == 0 and self.__Orientation.Z == True ):
            self.__infinte_length = True
        elif ( Nz <= 0 ):
            raise TypeError("the number of bin for z direction is not valid for cylinder filter with id " + str(filter_id) + ".")
        
        # get the low position 
        self.__low = Position(self.__origin.x() - pitch_x * 0.5, 
                              self.__origin.y() - pitch_y * 0.5, 
                              self.__origin.z())
        
        self.__x_index = 0
        self.__y_index = 0
        self.__z_index = 0
        if (self.__Real_Nx == 1):
            self.__x_index = 0
            self.__y_index -= 1
            self.__z_index -= 1
        
        if (self.__Real_Ny == 1):
            self.__x_index = 0
            self.__y_index = 0
            self.__z_index -= 1
        
        if (self.__Real_Nz == 1):
            self.__z_index = 0
    
    # function to reduce dimensions
    def __reduce_dimension(self, indices: list):
        if ( self.__Real_Nx == 1 and self.__Real_Ny == 1 and self.__Real_Nz == 1 ):
            return indices[0]
        
        if ( self.__infinte_length == True ):
            if ( self.__Nx == 1 and self.__Ny == 1 ):
                return indices[0]
        
        reduce_index = []
        if ( self.__Real_Nx > 1 ):
            reduce_index.append( indices[0] )
         
        if ( self.__Real_Ny > 1 ):
            reduce_index.append( indices[1] )
        
        if ( self.__Real_Nz > 1 ):
            reduce_index.append( indices[2] )
        
        return reduce_index
        
    def type_str(self):
        return "cylinder-filter"
    
    # get the indixes of the position
    def get_index(self, point: Position):
        mapped_r = self.__map_coordinate(point)
        
        nx = math.floor( (point.x() - self.__low.x()) / self.__pitch_x )
        ny = math.floor( (point.y() - self.__low.y()) / self.__pitch_y )
        if ( self.__infinte_length == True):
            nz = 0
        else:
            nz = math.floor( (point.z() - self.__low.z()) / self.__pitch_z )     
        
        
        new_origin_x = self.__origin.x() - self.__pitch_x * self.__Nx
        new_origin_y = self.__origin.y() - self.__pitch_y * self.__Ny
        indices = []
        if ( nx >=0 and nx < self.__Nx 
            and ny >= 0 and ny < self.__Ny
            and ((nz >= 0 and nz < self.__Nz) or self.__infinte_length ) ):
            
            if ( math.sqrt( (new_origin_x - mapped_r.x())**2 + (new_origin_y - mapped_r.y())**2 ) <= self.__radius + 1E-15 ):                    
                indices.append( nx )
                indices.append( ny )
                indices.append( nz )
            self.__map_index(indices)
            return (True, self.__reduce_dimension(indices))

        return (False, indices)

    def get_center(self, indices: list, is_map : bool = True)-> Position :
        new_origin_x = self.__origin.x()
        new_origin_y = self.__origin.y()
        new_origin_z = self.__origin.z()
        if ( self.__Real_Nx > 1 ):
            new_origin_x += self.__pitch_x * indices[self.__x_index]
        
        if ( self.__Real_Ny > 1 ):
            new_origin_y += self.__pitch_y * indices[self.__y_index]
        
        if( self.__Real_Nz > 1 ):
            new_origin_z += self.__pitch_z * indices[self.__z_index]
        
        if ( is_map == False ):
            return Position(new_origin_x, new_origin_y, new_origin_z)
        
        return self.__map_coordinate(new_origin_x, new_origin_y, new_origin_z)

    def get_scaled_radius_and_center(self, indices: list, r: Position)-> tuple :
        mapped_r = self.__map_coordinate(r)
        new_origin = self.get_center(indices, False)
        
        base = mapped_r.x() - new_origin.x()
        height = mapped_r.y() - new_origin.y()
        print(mapped_r)
        print(new_origin)
        scaled_r = math.sqrt(base* base + height * height) / self.__radius
        print(scaled_r)
        theta = math.atan2(height, base)
        if ( theta < 0. ):
            return (scaled_r, theta + 2*math.pi)
        
        return (scaled_r, theta)
        
    def z_min(self, indices: list) -> float :
        if ( self.__Real_Nz == 1 ):
            return self.__low.z()
        
        return self.__low.z() + self.__pitch_z * indices[self.__z_index]
    
    def z_max(self, indices: list) -> float :
        if ( self.__Real_Nz == 1 ):
            return self.__low.z() + self.__pitch_z
        
        return self.__low.z() + self.__pitch_z * (indices[self.__z_index] + 1)
    
    def dV(self):
        return math.pi * self.__radius * self.__radius * self.__pitch_z
        