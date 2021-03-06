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
        self.num_sources = len(body.panel_centre)
        self.strength = np.zeros(self.num_sources)
        self.coords = body.panel_centre
        self.coords[:, 1] = body.centreline  # Put sources on the centreline
        self.calc_sources()

        return None

    def calc_sources(self):
        """
        Calculate the strength of sources
        """
        # Calculate the strength of all the sources
        self.strength = self.calc_source_strength(self.body.mesh.normals,
                                                  [self.tank.U, 0, 0],
                                                  self.body.panel_area)

        # Remove sources above the waterline
        self.strength[np.argwhere(self.body.panel_centre[:, 2] > 0)] = 0

        # Remove sources with a negative y, relative to the centreline
        self.strength[np.argwhere(self.body.panel_centre[:, 1] -
                                  self.body.centreline > 0)] = 0

        return None

    def calc_source_strength(self, n, U, A):
        """
        Returns the strength of a source, given the normal vector, the onset
        free stream vector and the panel area
        """
        return (-1/(2 * np.pi))*np.dot(n, U) * A
