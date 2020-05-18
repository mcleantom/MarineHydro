# -*- coding: utf-8 -*-
"""
Created on Sun May 10 11:49:18 2020

@author: mclea
"""
import spam.hull as hull
import spam.sources as sources
import spam.tank as tank
import spam.wave_system as wave_system
import spam.stl_slicer as slicer
import numpy as np
# from mpl_toolkits import mplot3d
import matplotlib.pyplot as pyplot


def rw(x):
    t = tank.tank_properties()
    t.M = 300  # Maximum number of wave harmonics to calculate
    t.H = x  # The tank height
    t.B = 10  # The tank breadth
    t.U = (9.81/2)**0.5  # The flow speed
    t.reload()
    form_factor = 1.4
    h = hull.create_hull("data\\models\\5s.stl")
    h.mesh.translate([0, 0, 0])  # Translate X, Y, Z
    h.load_hull()
    s = sources.create_sources(h, t)
    waves = wave_system.surface_elevation(t, s)
#    print("The wave resistance = " + str(sum(waves.Rwm)) + " N")
    return sum(waves.Rwm)

speed_range = np.linspace(1,20, 100)
rrw = [rw(x) for x in speed_range]
pyplot.figure()
pyplot.plot(-speed_range, rrw)