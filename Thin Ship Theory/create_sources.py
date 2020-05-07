"""
Calculates the sources .

Classes:
    create_sources()

Methods:

Imports:
    numpy
"""
import numpy as np


class create_sources:
    """
    Creates all the sources for a hull

    Attributes:
        strength                    --  An array of source strengths
        coords                      --  The source's [X,Y,Z] coordinates

    Methods:
        calc_sources                -- Calculate the strength of the sources
    """

    def __init__(self, body, tank):
        """
        Initialise sources

        Inputs:
            body -- A create_hull object
            tank -- a tank object

        Outputs:
            A sources object
        """
        self.body = body
        self.tank = tank
        self.strength = np.zeros(len(body.panel_centre))
        self.coords = body.panel_centre
        self.coords[:, 0] = body.centreline  # Put sources on the centreline

        self.calc_sources()

    def calc_sources(self):
        """
        Calculate the strength of sources
        """
        self.strength = self.calc_source_strength(self.body.mesh.normals,
                                                  [self.tank.U, 0, 0],
                                                  self.body.panel_area)

        # ====================================================================
        #       Remove sources above waterline
        # ====================================================================
        for i in range(len(self.strength)):
            if self.body.panel_centre[i][2] > 0:
                self.strength[i] = 0

        # ====================================================================
        #       Remove sources with negative y (body assumed symetric,
        #       accounted for in the math)
        # ====================================================================
        for i in range(len(self.strength)):
            if self.body.panel_centre[i][1] < 0:
                    self.strength[i] = 0

#        print(self.strength)

    def calc_source_strength(self, n, U, A):
        """
        Returns the strength of a source, given the normal vector, the onset
        free stream vector and the panel area
        """
        return (-1/(2 * np.pi))*np.dot(n, U) * A
