"""
Creates a virtual hull.

Classes:
    create_hull()

Methods:

Imports:
    numpy, mesh from stl
"""

import numpy as np
from stl import mesh as stl_mesh


class create_hull:
    """
    Creates a hull object.

    Attributes:
        file_loc        --  File location of the hull
        mesh            --  stl object created from the file location
        panel_centre    --  An array of the coordinates of the centre of panels
        panel_area      --  An array of the areas of panels
        L               --  The length of the hull
        B               --  The breadth of the hull
        T               --  The draught of the hull
        D               --  The depth of the hull
        centreline      --  The Y coordinate of the centreline of the hull
        WSA             --  The wetted surface area of the hull

    Methods:
        load_hull       -- Reload the hull after a change in the mesh. (e.g. a translation)
    """

    def __init__(self, file):
        """
        Initialize a hull

        Inputs:
            file -- The location of an stl file

        Outputs:
            A hull object
        """

        self.file_loc = file
        self.mesh = stl_mesh.Mesh.from_file(self.file_loc)
        self.panel_centre = []
        self.panel_area = []
        self.L = 0  # Length
        self.B = 0  # Breadth
        self.T = 0  # Draft
        self.D = 0  # Depth
        self.centreline = 0  # Y coordinate of the centreline
        self.WSA = 0
        self.load_hull()

    def load_hull(self):
        """
        Reload the hull after a transformation of the mesh
        """
        self.calc_panel_centres()
        self.calc_panel_areas()
        self.calc_hull_parameters()

    def calc_panel_centres(self):
        """
        Calculate the centre of each panel
        """
        panel_x = (self.mesh.v0[:, 0]+self.mesh.v1[:, 0]+self.mesh.v2[:, 0])/3
        panel_y = (self.mesh.v0[:, 1]+self.mesh.v1[:, 1]+self.mesh.v2[:, 1])/3
        panel_z = (self.mesh.v0[:, 2]+self.mesh.v1[:, 2]+self.mesh.v2[:, 2])/3
        self.panel_centre = np.array([panel_x, panel_y, panel_z]).transpose()

    def calc_panel_areas(self):
        """
        Calculate the area of each panel
        """
        AB = self.mesh.v1-self.mesh.v0
        AC = self.mesh.v2-self.mesh.v0
        self.panel_area = np.sum(0.5*abs(np.cross(AB, AC)), axis=1)

    def calc_hull_parameters(self):
        """
        Calculate the hull length, breadth, draught, depth and the centreline
        """
        self.L = np.amax(self.mesh.x) - np.amin(self.mesh.x)
        self.B = np.amax(self.mesh.y) - np.amin(self.mesh.y)
        self.T = abs(np.amin(self.mesh.z))
        self.D = np.amax(self.mesh.z) - np.amin(self.mesh.z)
        self.centreline = np.average(self.mesh.y)
        self.WSA = sum(self.panel_area * (self.panel_centre[:, 2] < 0))
