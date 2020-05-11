"""
Classes:
    wave_resistance()

Methods:

Imports:
    numpy
"""
import numpy as np
from scipy import optimize

class surface_elevation:
    """
    Creates a wave system for a given tank and distribution of sources within
    the tank

    Attributes:
        m       -- Wave harmonics
        etam    -- Wave amplitude coefficient
        num     -- Wave amplitude coefficient
        km      -- Wave number
        thetam  -- Wave angle
        Rwm     -- Wave reistance

    Methods:
        calc_surface_elevation  -- Calculates the surface elevations etam, num
        calc_Rwm                -- Calculates the wave resistance due to component m
    """
    def __init__(self, tank, sources):
        """
        Initialises the wave_system class, calculates the wave components and
        wave reistances for harmonics m.
        """
        self.tank = tank
        self.sources = sources
        self.etam = self.blank_list()  # Wave amplitude coefficieints = 2*An*cos(en)
        self.num = self.blank_list()  # Wave amplitude coefficient = 2*An*sin(en)
        self.km = self.blank_list()
        self.thetam = self.blank_list()
        self.Rwm = self.blank_list()
        self.m = np.arange(0, tank.M)
        self.wave_components()
        self.calc_surface_elevation()
        self.calc_Rwm()
        return None

    def calc_Rwm(self):
        """
        Calculates the wave reistance values at each component
        """
        coeff = ((16*np.pi*self.tank.U)/(self.tank.B*self.tank.g))
        self.zetam = self.etam*self.num
        # m = 0
        bracket_term = (1 - ((2*self.km[0]*self.tank.H) /
                             (np.sinh(2*self.km[0]*self.tank.H))))
        self.Rwm[0] = coeff*self.zetam[0]*bracket_term

        # m >= 1
        self.Rwm[1::] = self.calc_rw_summation()

        return None

    def calc_rw_summation(self):
        """
        Calculates the summation component in Rw (m>=1)
        """
        bracket_term = 1 + ((2*self.km[1::]*self.tank.H) /
                            (np.sinh(2*self.km[1::]*self.tank.H)))
        cos_term = (np.cos(self.thetam[1::])**2)/2
        return self.zetam[1::]*(1-cos_term*bracket_term)

    def blank_list(self):
        """
        Returns a np array with zeros of a length number of wave components
        """
        return np.zeros(self.tank.M)

    def calc_surface_elevation(self):
        """
        Calculate the surface elevation coefficients over a range of wave
        components for a given distribution of sources in a tank
        """

        coeff = ((16*np.pi*self.tank.U)/(self.tank.B*self.tank.g))

        top_frac = ((self.km[0] + self.km*np.cos(self.thetam)**2))
        sechKH = (1/(np.cosh(self.tank.H*self.km)))
        bottom_frac = (1 + np.sin(self.thetam)**2 - self.km[0]*self.tank.H*sechKH**2)
        frac_term = top_frac/bottom_frac
        summation_term = self.calc_elevation_summation()
        wave_comps = coeff*frac_term*summation_term
        wave_comps[0] = wave_comps[0] * 0.5
        self.etam = wave_comps[0]
        self.num = wave_comps[1]

        return None

    def calc_elevation_summation(self):
        """
        Calculates the summation component in the formula for etam and num
        """
        x_sigma, y_sigma, z_sigma = self.separate_coords(self.sources.coords)
        m_matrix = np.array([self.m, ] * len(y_sigma)).T
        k_matrix = np.array([self.km, ] * len(y_sigma)).T
        theta_matrix = np.array([self.thetam, ] * len(y_sigma)).T

        sigma_term = self.sources.strength
        exp_term = np.exp(-1*k_matrix*self.tank.H)
        cosh_term = np.cosh((k_matrix*(self.tank.H+z_sigma)))
        matrix_term = np.array([np.cos(k_matrix*x_sigma*np.cos(theta_matrix)),
                                np.sin(k_matrix*x_sigma*np.cos(theta_matrix))])
        mpiyoB = m_matrix*np.pi*y_sigma/self.tank.B
        cos_sin_term = np.zeros(mpiyoB.shape)
        cos_sin_term[::2] = np.cos(mpiyoB[::2])
        cos_sin_term[1::2] = np.sin(mpiyoB[1::2])

        summation_terms = (sigma_term*exp_term*cosh_term*matrix_term*cos_sin_term)

        etan = np.sum(summation_terms[0], axis=1)
        nun = np.sum(summation_terms[1], axis=1)

        return np.array([etan, nun])

    def separate_coords(self, coords):
        """
        Returns separated coordinates for a given list of coordinates
        """
        return coords[:, 0].T, coords[:, 1].T, coords[:, 2].T

    def wave_components(self):
        """
        Calculates the wave components of for a range of wave components from 0
        to a maximum value M, using a newton raphson method

        Inputs:
            M   -- Maximum wave harmonic number
        """
        X = self.m * np.pi/self.tank.B  # Wall condition
        kguess = self.tank.k0 * (1 + (1 + (2 * X/self.tank.k0))**0.5) / 2  # Initial guess
        self.km = np.array([optimize.newton(self.f, k, args=(m,)) for k, m in zip(kguess, self.m)]) # Km = The roots of the equation f(x)
        self.thetam = np.arcsin(((2*np.pi*self.m)/self.tank.B)/self.km)

        return None

    def f(self, x, n,):
        """
        Returns f(x) = x**2-k0*yn*tanh(yn*H)+((2*m*pi)/B)^2

        Inputs:
            x   --  The wave number yn
            n   --  The wave harmonic to calculate for
        """
        return x**2 - self.tank.k0*x*np.tanh(x*self.tank.H) - 1*((2*n*np.pi)/(self.tank.B))**2
