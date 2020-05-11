# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:20:50 2020

@author: mclea
"""
import numpy as np

class resistance():
    """
    """
    def __init__(self, wave_system):
        """
        """
        self.Rw = []
        self.Rp = []
        self.tank = wave_system.tank
#        self.boat = wave_system.hull
        return None

    def calc_wave_resistance(self):
        """
        """
        coeff = (16*np.pi*self.tank.U)/(self.tank.B*self.tank.g)
        z0 = self.wave_system.etan[0]**2 + self.wave_system.nun[0]**2
        print(coeff*z0)
        return None

    def calc_summation(self):
        """
        """
        return None