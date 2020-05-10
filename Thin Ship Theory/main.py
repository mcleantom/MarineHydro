# -*- coding: utf-8 -*-
"""
Created on Tue May  5 14:26:52 2020

@author: mclea
"""
import spam.hull as hull
import spam.resistance as resistance
import spam.sources as sources
import spam.tank as tank
import spam.wave_system as wave_system
import numpy as np

tank = tank.tank_properties()
tank.M = 5  # Maximum number of wave harmonics to calculate
tank.H = 2  # The tank height
tank.B = 10  # The tank breadth

hull = hull.create_hull("data\\models\\5s.stl")
sources = sources.create_sources(hull, tank)

waves = wave_system.surface_elevation(tank, sources)
waves.wave_components(10)
print(waves.km)
print(np.degrees(waves.thetam))
