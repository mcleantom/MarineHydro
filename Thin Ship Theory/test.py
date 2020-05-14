# -*- coding: utf-8 -*-
"""
Created on Sun May 10 11:49:18 2020

@author: mclea
"""
import spam.hull as hull
import spam.sources as sources
import spam.tank as tank
import spam.wave_system as wave_system

tank = tank.tank_properties()
tank.M = 200  # Maximum number of wave harmonics to calculate
tank.H = 1  # The tank height
tank.B = 20  # The tank breadth
tank.U = -10  # The flow speed

hull = hull.create_hull("data\\models\\5s.stl")
hull.mesh.translate([0, 0, 0])
hull.load_hull()
sources = sources.create_sources(hull, tank)

waves = wave_system.surface_elevation(tank, sources)
print("The wave resistance = " + str(sum(waves.Rwm)) + " N")

waves.plot()
