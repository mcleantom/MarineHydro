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
tank.M = 150  # Maximum number of wave harmonics to calculate
tank.H = 2  # The tank height
tank.B = 10  # The tank breadth
tank.U = 1  # The flow speed

hull = hull.create_hull("data\\models\\cube.stl.txt")
sources = sources.create_sources(hull, tank)
test_source = sources
test_source.num_sources = 1
test_source.strength = np.array([1])
test_source.coords = np.array([[2, 2, -1]])

waves = wave_system.surface_elevation(tank, sources)
print("The wave resistance = " + str(sum(waves.Rwm)) + " N")

nx, ny = (200, 100)
num_points = nx*ny
x = np.linspace(-2, 10, nx)
y = np.linspace(-5, 5, ny)
xim = waves.xim
etam = waves.etam
km = waves.km
thetam = waves.thetam
xx, yy = np.meshgrid(x, y)
xx = np.reshape(xx.flatten(), (num_points,1))
yy = yy.flatten()
yy = np.reshape(yy.flatten(), (num_points,1))
m = waves.m

term1 = xim*np.cos(xx*km*np.cos(thetam))
term2 = etam*np.sin(xx*km*np.cos(thetam))
term3 = np.cos((m*np.pi*yy)/tank.B)
z_x_y = (term1*term2)*term3
z_x_y = np.sum(z_x_y, axis=1)
z_x_y[np.argwhere(xx < 0)] = 0
z_x_y = z_x_y.reshape(ny, -1)
xx = xx.reshape(ny, -1)
yy = yy.reshape(ny, -1)

fig = plt.figure()
h = plt.contourf(xx, yy, z_x_y)
plt.show()