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
tank.M = 5  # Maximum number of wave harmonics to calculate
tank.H = 2  # The tank height
tank.B = 10  # The tank breadth

hull = hull.create_hull("data\\models\\5s.stl")

sources = sources.create_sources(hull, tank)
Rw = resistance.wave_resistance(sources, tank)

for i, value in enumerate(Rw.RWm):
    print("Rw for harmonic " + str(i) + " is " + str(value) + " N")

print(sum(Rw.RWm))
