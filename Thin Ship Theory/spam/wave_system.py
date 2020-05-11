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
    
    """
    def __init__(self, tank, sources):
        """
        """
        self.tank = tank
        self.sources = sources
        self.etan = self.blank_list()  # Wave amplitude coefficieints = 2*An*cos(en)
        self.nun = self.blank_list()  # Wave amplitude coefficient = 2*An*sin(en)
        self.km = self.blank_list()
        self.thetam = self.blank_list()
        self.m = np.arange(0,tank.M)
#        print(self.etan)
        self.wave_components()
        self.calc_surface_elevation()

        return None

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
        summation_term = self.calc_summation()
#
#        wave_comps = coeff*frac_term*summation_term
#        wave_comps[0] = wave_comps[0] * 0.5
#        self.etan = wave_comps[0]
#        self.nun = wave_comps[1]
        
        return None

    def calc_summation(self):
        """
        
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
        
        self.etan = np.sum(summation_terms[0], axis=1)
        self.nun = np.sum(summation_terms[1], axis=1)
        
        return [self.etan, self.nun]


    def separate_coords(self, coords):
        return coords[:,0].T, coords[:,1].T, coords[:,2].T

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
