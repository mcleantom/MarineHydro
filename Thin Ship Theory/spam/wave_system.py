"""
Classes:
    surface_elevation()

Methods:

Imports:
    numpy, optimize from scipy
"""
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

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
        self.xim = self.blank_list()  # Wave amplitude coefficieints = 2*An*cos(en)
        self.etam = self.blank_list()  # Wave amplitude coefficient = 2*An*sin(en)
        self.km = self.blank_list()  # Wave numers
        self.thetam = self.blank_list()  # Wave angles
        self.Rwm = self.blank_list()  # Components of wave resistance
        self.m = np.arange(0, tank.M)  # Array of wave components
        self.wave_components()
        self.calc_surface_elevation()
        self.calc_Rwm()
        return None

    def calc_Rwm(self):
        """
        Calculates the wave reistance values at each component
        """
        coeff = ((16*np.pi*self.tank.U)/(self.tank.B*self.tank.g))
        self.zetam = self.xim*self.etam
        # m = 0
        bracket_term = (1 - ((2*self.km[0]*self.tank.H) /
                             (np.sinh(2*self.km[0]*self.tank.H))))
        self.Rwm[0] = coeff*self.zetam[0]*bracket_term

        # m >= 1
        self.Rwm[1::] = coeff*self.calc_rw_summation()

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
#        print(frac_term)
        summation_term = self.calc_elevation_summation()
        wave_comps = coeff*frac_term*summation_term
        wave_comps[0] = wave_comps[0] * 0.5
        self.xim = wave_comps[0]
        self.xim = np.nan_to_num(self.xim)
        self.etam = wave_comps[1]
        self.etam = np.nan_to_num(self.etam)
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

        # Make a matrix of the cos(mpiy/B) and sin(mpiy/B) terms
        matrix_term = np.array([np.cos(k_matrix*x_sigma*np.cos(theta_matrix)),
                                np.sin(k_matrix*x_sigma*np.cos(theta_matrix))])
        mpiyoB = m_matrix*np.pi*y_sigma/self.tank.B
        cos_sin_term = np.zeros(mpiyoB.shape)
        cos_sin_term[::2] = np.cos(mpiyoB[::2])
        cos_sin_term[1::2] = np.sin(mpiyoB[1::2])

        summation_terms = (sigma_term*exp_term*cosh_term*matrix_term*cos_sin_term)

        xi = np.sum(summation_terms[0], axis=1)
        eta  = np.sum(summation_terms[1], axis=1)

        return np.array([xi, eta])

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
        return (x**2 - self.tank.k0*x*np.tanh(x*self.tank.H) -
                1*((2*n*np.pi)/(self.tank.B))**2)

    def calc_wave_height(self, nx=500, minx=-1, maxx=15,
                         ny=100, miny=-3, maxy=3):
        """
        Plot the wave profile
        Inputs:
            nx      --  Number of points in the x axis
            minx    --  Minumum value of x
            maxx    --  Maximum value of x
            ny      --  Number of points in the y axis
            miny    --  Minimum value of y
            maxy    --  Maximum value of y
        """
        # Set the plot of the tank
        miny = -1*self.tank.B/2
        maxy = -1*miny
        num_points = nx*ny
        x = np.linspace(minx, maxx, nx)
        y = np.linspace(miny, maxy, ny)
        self.xx, self.yy = np.meshgrid(x, y)
        self.xx = np.reshape(self.xx.flatten(), (num_points, 1))
        self.yy = self.yy.flatten()
        self.yy = np.reshape(self.yy.flatten(), (num_points, 1))

        # Calculate the wave height
        term1 = self.xim*np.cos(self.xx*self.km*np.cos(self.thetam))
        term2 = self.etam*np.sin(self.xx*self.km*np.cos(self.thetam))
        term3 = np.cos((self.m*np.pi*self.yy)/self.tank.B)
        self.z_x_y = (term1*term2)*term3
        self.z_x_y = np.sum(self.z_x_y, axis=1)
        self.z_x_y[np.argwhere(self.xx < 0)] = 0
        self.z_x_y = self.z_x_y.reshape(ny, -1)
        self.xx = self.xx.reshape(ny, -1)
        self.yy = self.yy.reshape(ny, -1)

    def plot(self, x, y):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        cf = ax.contourf(self.xx, self.yy, self.z_x_y)
        plt.colorbar(cf, ax=ax)
        plt.plot(x, y, 'r.')
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("Wave pattern for a speed of " + str(self.tank.U) + " m/s")
