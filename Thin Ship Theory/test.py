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
# import matplotlib.pyplot as pyplot

tank = tank.tank_properties()
tank.M = 300  # Maximum number of wave harmonics to calculate
tank.H = 10  # The tank height
tank.B = 10  # The tank breadth
tank.U = 4  #(9.81/2)**0.5  # The flow speed
tank.reload()
hull = hull.create_hull("data\\models\\cube.stl.txt")
hull.mesh.translate([0, 0, 0.0])  # Translate X, Y, Z
hull.load_hull()
sources = sources.create_sources(hull, tank)

print(hull.WSA)
#waves = wave_system.surface_elevation(tank, sources)
#print("The wave resistance = " + str(sum(waves.Rwm)) + " N")
#
#waves.calc_wave_height()
#plane_normal = [0, 0, 1]  # The plane normal vector
#plane_origin = [0, 0, 0]  # A point on the plane
#plane = [plane_origin, plane_normal]  # A plane facing right
#
#hull_slice = slicer.make_slice(hull.mesh, plane)
#
#x = np.flip(hull_slice.slice_points[:, 0])
#y = np.flip(hull_slice.slice_points[:, 1])
#waves.plot(x, y)

#print(waves.xim)

#plane_normal = [0, 0, 1]  # The plane normal vector
#plane_origin = [0, 0, 0]  # A point on the plane
#plane = [plane_origin, plane_normal]  # A plane facing right
#
#hull_slice = slicer.make_slice(hull.mesh, plane)
#
#x = np.flip(hull_slice.slice_points[:, 0])
#y = np.flip(hull_slice.slice_points[:, 1])
#
## Create a new plot
#figure = pyplot.figure()
#axes = mplot3d.Axes3D(figure)

# Load the STL files and add the vectors to the plot
#axes.add_collection3d(mplot3d.art3d.Poly3DCollection(hull.mesh.vectors))
#axes.contour(waves.xx, waves.yy, waves.z_x_y, cmap="RdYlBu")
# Auto scale to the mesh size
#scale = hull.mesh.points.flatten(-1)
#axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
#pyplot.show()
