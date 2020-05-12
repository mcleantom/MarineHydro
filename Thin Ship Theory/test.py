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
import matplotlib.pyplot as plt
import scipy.ndimage
#from mpl_toolkits import mplot3d

tank = tank.tank_properties()
tank.M = 100  # Maximum number of wave harmonics to calculate
tank.H = 2  # The tank height
tank.B = 10  # The tank breadth
tank.U = 15  # The flow speed

hull = hull.create_hull("data\\models\\5s.stl")
sources = sources.create_sources(hull, tank)
test_source = sources
test_source.num_sources = 1
test_source.strength = np.array([1])
test_source.coords = np.array([[2, 2, -1]])

waves = wave_system.surface_elevation(tank, sources)
print("The wave resistance = " + str(sum(waves.Rwm)) + " N")

nx, ny = (300,300)
x = np.linspace(-2,10,nx)
y = np.linspace(-5,5,ny)
xim = waves.xim
etam = waves.etam
km = waves.km
thetam = waves.thetam
xx, yy = np.meshgrid(x,y)
z_x_y = np.zeros((nx, ny))
m = waves.m

@np.vectorize
def return_elevation(xx, yy):
    summation = 0
    if xx>0:
        for i in m:
            term1 = xim[i] * np.cos(xx * km[i] * np.cos(thetam[i]))
            term2 = etam[i] * np.sin(xx * km[i] * np.cos(thetam[i]))
            term3 = np.cos((i*np.pi*yy)/(tank.B))
            summation += (term1+term2)*term3
    return summation

z_x_y = return_elevation(xx,yy)

fig = plt.figure()
h = plt.contourf(xx, yy, z_x_y)
plt.show()