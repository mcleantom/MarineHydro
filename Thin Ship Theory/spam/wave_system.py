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
        summation_term = self.calc_summation().T[0]

        wave_comps = coeff*frac_term*summation_term
        wave_comps[0] = wave_comps[0] * 0.5
        self.etan = wave_comps[0]
        self.nun = wave_comps[1]
        print(wave_comps)

    def calc_summation(self):
        """
        
        """
        return_values = np.zeros((len(self.m),2,1))
        
        for i in self.m:
            summation = 0
            for j in range(self.sources.num_sources):
                source_strength = self.sources.strength[j]
                x,y,z = self.separate_coords(self.sources.coords[j])
                exp_term = np.exp(self.km[i]*z)
                cosh_term = np.cosh(self.km[i]*(self.tank.H+z))
                matrix_term = np.array([[np.cos(self.km[i] * x *
                                            np.cos(self.thetam[i]))],
                                        [np.sin(self.km[i] * x *
                                            np.sin(self.thetam[i]))]])
                if i % 2 == 0:
                    multiplier = np.cos((i*np.pi*y)/self.tank.B)
                else:
                    multiplier = np.sin((i*np.pi*y)/self.tank.B)
                summation += source_strength*exp_term*cosh_term*matrix_term*multiplier
                
            return_values[i] = summation

        return return_values


    def separate_coords(self, coords):
        return coords[0], coords[1], coords[2]

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
