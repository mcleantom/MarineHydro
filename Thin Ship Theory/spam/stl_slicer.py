# -*- coding: utf-8 -*-
"""
This class takes a stl object and can slice it at any plane to create a 2D
cross section.

Classes:
    make_slice()

Methods:

Imports:
    numpy, mesh from stl
"""

import numpy as np


class make_slice:
    """
    Slice class

    Attributes:
        slice_points    --  An array of all the points in a slice of the stl
                            object
    """

    def __init__(self, obj, plane):
        """Initialise a slice of an stl object at a plane

        Inputs:
            obj     --  An stl.mesh object
            plane   --  A plane array in the format [P,N] where P is a position
                        on the plane and N is the normal vector
        Outputs:
            A slice of a hull

        Example:
            hull_slice = make_slice(stl_object, [[0,0,0],[1,0,0]]) # make a slice of the stl file
        """
        self.obj = obj
        self.plane = plane
        self.plane_origin = plane[0]
        self.plane_normal = plane[1]
        self.distance_from_plane = np.dot((self.obj.vectors-self.plane[0]),
                                          self.plane[1])
        self.slice_points = []

        self.calculate_points()
        
        return None

    def calculate_points(self):
        """Calculate the points on a slice of an stl file
        """

        v0 = np.vstack(np.sign(self.distance_from_plane[:, 0] *
                               self.distance_from_plane[:, 1]))
        v1 = np.vstack(np.sign(self.distance_from_plane[:, 1] *
                               self.distance_from_plane[:, 2]))
        v2 = np.vstack(np.sign(self.distance_from_plane[:, 2] *
                               self.distance_from_plane[:, 0]))

        goes_through = np.concatenate([v0, v1, v2], 1)
        self.num_points = (np.sum(goes_through == 0) * 2 +
                           np.sum(goes_through == -1))

        self.slice_points = np.zeros((self.num_points, 3))
        curr_point = 0

        for i in range(len(goes_through)):
            for x in range(len(goes_through[i])):
                if goes_through[i][x] == -1:
                    # The line goes through the plane
                    # There is a point which lies in the plane, which is r
                    # amounts of vector V away from point P0
                    p0 = self.obj.vectors[i][x]
                    if x == 0 or x == 1:
                        V = self.obj.vectors[i][x]-self.obj.vectors[i][x+1]
                    else:
                        V = self.obj.vectors[i][x]-self.obj.vectors[i][0]
                    d = (np.dot((self.plane_origin-p0), self.plane_normal) /
                         (np.dot(V, self.plane_normal)))
                    P = p0 + d*V
                    self.slice_points[curr_point] = P
                    curr_point += 1

                if goes_through[i][x] == 0:
                    if x == 0 or x == 1:
                        self.slice_points[curr_point] = self.obj.vectors[i][x]
                        self.slice_points[curr_point+1] = self.obj.vectors[i][x+1]
                    else:
                        self.slice_points[curr_point] = self.obj.vectors[i][x]
                        self.slice_points[curr_point+1] = self.obj.vectors[i][0]

        self.slice_points = np.round(self.slice_points, decimals=5)
        self.slice_points = np.unique(self.slice_points, axis=0)
        
        return None
        
    