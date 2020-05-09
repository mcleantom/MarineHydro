# -*- coding: utf-8 -*-
"""
Created on Tue May  5 14:26:52 2020

@author: mclea
"""

import spam.hull as hull
import spam.resistance as resistance
import spam.sources as sources
import spam.tank as tank

tank = tank.tank_properties()
tank.M = 5
tank.H = 2
tank.B = 10

hull = hull.create_hull('5s.stl')
sources = sources.create_sources(hull, tank)
Rw = resistance.wave_resistance(sources, tank)

for i, value in enumerate(Rw.RWm):
    print("Rw for harmonic " + str(i) + " is " + str(value) + " N")

print(sum(Rw.RWm))
