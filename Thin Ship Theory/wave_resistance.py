"""
Classes:
    wave_resistance()

Methods:

Imports:
    numpy, optimize from scipy
"""
import numpy as np
from scipy import optimize


class wave_resistance:
    """
    Creates a class which deals with calculating wave resistances

    Attributes:
        RWm         --  An array of Rw values at each wave harmonic
        elevation   --  Wave elevation terms
    """

    def __init__(self, sources, tank):
        """
        Initialise sources

        Inputs:
            sources -- A sources object
            tank -- a tank object

        Outputs:
            A wave-resistance class
        """
        self.sources = sources
        self.tank = tank
        self.elevation = []
        self.RWm = np.zeros(tank.M)
        self.kBar = 0

        self.calc_Rwm()

    def calc_Rwm(self):
        """
        Calculate the wave resistance values
        """
        self.kbar = self.tank.g/(self.tank.U**2)

        coeff = 0.25*self.tank.rho*self.tank.g*self.tank.B

        k0, theta0 = self.wave_components(0)
        eta0, nu0 = self.elevation_terms(0)
        zeta0_squared = eta0**2 + nu0**2

        zeroth_component_to_rw = self.zeroth_wave_component(k0, zeta0_squared)
        self.RWm[0] = coeff*zeroth_component_to_rw

        for i in np.arange(1, self.tank.M):
            km, thetam = self.wave_components(i)
            etam, num = self.elevation_terms(i)
            zetam_squared = etam**2 + num**2

            mth_component_to_rw = self.mth_wave_component(km, thetam, zetam_squared)
            self.RWm[i] = coeff*mth_component_to_rw

    def zeroth_wave_component(self, k0, zeta0_squared):
        """
        Calculates the component to Rw for the 0th wave component
        """
        if (2*k0*self.tank.H) > 50:
            frac_term = 0
        else:
            top_frac = 2*k0*self.tank.H
            bottom_frac = np.sinh(2*k0*self.tank.H)
            frac_term = top_frac/bottom_frac

        return zeta0_squared*(1-frac_term)
    
    def mth_wave_component(self, km, thetam, zetam_squared):
        """
        Calculates the component to Rw for the mth wave component
        """
        frac1 = (np.cos(thetam)**2)/2

        if (2*km*self.tank.H) > 10:
            frac2 = 0
        else:
            top_frac = (2*km*self.tank.H)
            bottom_frac = (np.sinh(2*km*self.tank.H))
            frac2 = top_frac/bottom_frac

        return zetam_squared*(1-(frac1*(1+frac2)))
# =============================================================================
# POSSIBLE CALCULATING WAVE COMPONENTS IS WRONG. READ PAPER AND FIX.
# =============================================================================

    def wave_components(self, m):
        """
        Calculates the wave components of a given harmonic m using a newton
        raphson method of root finding.

        Inputs:
            m   --  Wave harmonic number
        """
        X = m * np.pi/self.tank.B  # Wall condition
        k = self.tank.k0 * (1 + (1 + (2 * X/self.tank.k0))**0.5) / 2  # Initial guess
        km = optimize.newton(self.f, k, args=(m,))
        thetan = np.arcsin(((2*np.pi*m)/self.tank.B)/km)

        return [km, thetan]

    def f(self, x, n,):
        """
        Returns f(x) = x**2-k0*yn*tanh(yn*H)+((2*m*pi)/B)^2

        Inputs:
            x   --  The wave number yn
            n   --  The wave harmonic to calculate for
        """
        return x**2 - self.tank.k0*x*np.tanh(x*self.tank.H) - 1*((2*n*np.pi)/(self.tank.B))**2

    def elevation_terms(self, m):
        """
        Calculates the wave elevation terms for a given harmonic m.

        Inputs:
            m   --  Wave harmonic number
        """
        coeff = (16*np.pi*self.tank.U)/(self.tank.B*self.tank.g)

        km, thetam = self.wave_components(m)

        first_frac = self.calc_first_frac(km, thetam)
        summation = self.calc_summation_term(m, km, thetam)

        if m == 0:              # If the 0th wave component
            coeff = coeff*0.5   # Half the output

        return coeff*first_frac*summation

    def calc_first_frac(self, km, thetam):
        """
        Calculates the first fraction term in the wave elevation equation
        returns ((Kbar + km*cos^(thetam))/(1+sin^2(thetam)-Kbar*H*sech^2(km*H)))
        """
        top_frac = (self.kbar+km*(np.cos(thetam)**2))

        if km*self.tank.H > 20:
            sechKH = 0
        else:
            sechKH = 1/np.cosh(km*self.tank.H)

        bottom_frac = (1 + (np.sin(thetam)**2) - self.kbar*self.tank.H*(sechKH**2))
        first_frac = top_frac/bottom_frac
        return first_frac

    def calc_summation_term(self, m, km, thetam):
        """
        Calculates the summation term of the wave elevation equation
        """
        summation = 0

        for i in range(len(self.sources.strength)):
            strength_i = self.sources.strength[i]
            exp_term = np.exp(-1*km * self.tank.H)
            cosh_term = np.cosh(km*(self.tank.H * self.sources.coords[i][2]))
            matrix_term = np.array([[np.cos(km * self.sources.coords[i][0] *
                                            np.cos(thetam))],
                                    [np.sin(km*self.sources.coords[i][0] *
                                            np.cos(thetam))]])
            if m % 2 == 0:  # Even m
                multiplier = np.cos((m*np.pi*self.sources.coords[i][1]) /
                                    self.tank.B)
            else:  # Odd m
                multiplier = np.sin((m*np.pi*self.sources.coords[i][1]) /
                                    self.tank.B)

            summation += strength_i*exp_term*cosh_term*matrix_term*multiplier

        return summation
