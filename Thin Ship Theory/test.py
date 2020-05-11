# -*- coding: utf-8 -*-
"""
Created on Sun May 10 11:49:18 2020

@author: mclea
"""
import spam.hull as hull
import spam.sources as sources
import spam.tank as tank
import spam.wave_system as wave_system
import numpy as np

tank = tank.tank_properties()
tank.M = 5  # Maximum number of wave harmonics to calculate
tank.H = 2  # The tank height
tank.B = 5  # The tank breadth

hull = hull.create_hull("data\\models\\5s.stl")
sources = sources.create_sources(hull, tank)
#test_source = sources
#test_source.num_sources = 1
#test_source.strength = np.array([0.5])
#test_source.coords = np.array([[2, 2, -1]])

waves = wave_system.surface_elevation(tank, sources)

#print(waves.etan)
#print(waves.nun)