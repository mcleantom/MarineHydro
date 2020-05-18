"""
Creates a virtual tank.

Classes:
    tank_properties

Methods:

Imports:

"""


class tank_properties:
    """
    Creates a tank

    Attributes:
        g -- Acceleration due to gravity
        c -- Ship speed
        K0 -- Fundimental wave number (g/c^2)
        U -- Ship speed
        B -- Tank width
        H -- Tank height
        M -- Maximum number of wave harmonics
        rho -- Water density
        mu -- Water viscosity
    """

    def __init__(self):
        """ Initialises a tank

        Inputs:
            None

        Outputs:
            A tank
        """
        self.g = 9.81  # Acceleration due to gravity
#        self.c = (self.g/2)**0.5  # Ship speed
#        self.k0 = self.g/self.c**2
        self.U = 2  # Water speed (m/s)
        self.B = 10  # Tank width (m)
        self.H = 20  # Tank depth
        self.M = 2  # Max wave harmonic
        self.rho = 1000  # Water density
        self.mu = 1.139e-06  # viscosity
        self.reload()
        return None
    
    def reload(self):
#        self.c = -1*self.U
        self.k0 = self.g/self.U**2