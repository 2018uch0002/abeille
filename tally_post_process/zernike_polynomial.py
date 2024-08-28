#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 19:40:30 2024

@author: parthsingh
"""

import math
import numpy

# =============================================================================
# Zernike Polynomial 
# =============================================================================
class ZernikePolynomials:
    '''
    ZernikePolynomials class take order as a parameter and evaluate all the coefficients up to given order.
    '''
    def __init__(self, order: int):

        if ( isinstance(order, int) or isinstance(order, numpy.uint64)):
            self.__order = order
        else:
            raise ValueError("Zernike Polynomial expect the order to be an int.")
        
        # to store the coefficients
        self.__Zr_coefficents = []
        # to store the indication marks for odd polynomial (1) or even even polynomial(0)
        self.__Zr_type = [] 
        # store the n
        self.__n = []
        # store the m
        self.__m = []
        
        # evaluate and store all the coefficients
        for i in range(0, order+1):
            (n, l) = self.__get_n_and_l(i)
            m = abs(l)
            # store the odd or even indication
            if ( l < 0 ):
                self.__Zr_type.append(1)
            elif ( l >= 0 ) :
                  self.__Zr_type.append(0)
            
            # store the n and m
            self.__m.append(m)
            self.__n.append(n)
            
            # max power 
            max_k = int( (n-m)/2 )
            sign_change = 1
            if ( max_k % 2 == 1 ):
                sign_change = -1
            
            for j in range(0, max_k+1):
                k = max_k - j
                avg_n_m = int( (n+m) / 2 )
                mid_n_m = int( (n-m) / 2 )
                numerator = self.__factorial( n - k )
                denominator = self.__factorial( k ) * self.__factorial( avg_n_m - k ) * self.__factorial( mid_n_m - k )
                coeff = sign_change * numerator / denominator
                
                self.__Zr_coefficents.append( coeff )
                sign_change *= -1
            
            # store the 'n' which will correponds to highest given order
            self.__max_n = n
            
    # get the order
    def order(self):
        '''
        return the highest possible order in the class instance

        '''
        return self.__order
    
    # evaluate all the orders
    def eval_zernike(self, r : float , theta : float):
        '''
        provides the list of values, evaluated zernike at ach order.
        '''
        
        # check for correct scaled_r and theta
        if ( r > 1 ):
            raise ValueError( "the given \"scaled_r\" must be between [0, 1].")
        if ( abs(theta) > math.pi ):
            raise ValueError("the given \"theta\" must be between [-pi, pi] .")
        
        # first evaluate all the powers
        r_powers = [1.]
        r_k = 1.
        for i in range(0, self.__max_n + 1):
            r_k *= r
            r_powers.append( r_k )
        
        # evaluate each order polynomial value and store into the list
        it_coeff = 0
        zr_values = []
        for i in range(0, self.__order + 1):
            m = self.__m[i]
            n = self.__n[i]
            order_coeff_size = int( (n-m)/2 + 1 )
            zr_type = self.__Zr_type[i]
            value = 0.
            k = m
            for j in range(0, order_coeff_size):
                value += self.__Zr_coefficents[it_coeff] * r_powers[k]
                k += 2
                it_coeff += 1
            if ( zr_type == 0 ):
                value *= math.cos( m * theta )
            elif ( zr_type == 1 ):
                value *= math.sin( m * theta )
            
            zr_values.append(value)
        
        return zr_values
        
    
    # orthonormalisation constant defined as inverse of norm-2
    def orthonormalisation_constant(self, order: int):
        '''
        return the orthonormalisation constant defined as inverse of norm-2
        note that, it does not include \"pi\" in the multipication."
        '''
        (n, l) = self.__get_n_and_l(order)
        if ( l == 0 ):
            return n+ 1.
        else:
            return 2. * ( n + 1. )
        
    
    
    # the method to get unique int pair of 'n' and 'l' 
    # order = 0.5 * [n*(n+2) + l]
    def __get_n_and_l(self, order):
        n = 0
        while( n <= order ):
            l = 2 * order - n * (n+2)
            if ( abs(l) <= n ):
                return (n, l)
            
            n += 1
        raise ValueError("A unique pair of 'n' and 'l'  is not found for order  = " + str (order) )
        return (0, 0)
    
    # factorial function 
    def __factorial(self, N):
        if (N == 0 or N == 1):
            return 1.
        else:
            
            return N * self.__factorial(N-1)
      
        
        