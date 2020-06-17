#!/usr/bin/env python3
import csv                               ## for reading in our data files
import logging                           ## for orderly print output
import numpy as np                       ## for handling of data
import sys                               ## useful system calls (used to exit cleanly)
import os                                ## path manipulations
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ipywidgets import interact, interactive, fixed, widgets


def GaussFunc(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def LineFunc(x, k, m):
    return k*x+m


class Gauss:
    """A class to hold coefficients for Gaussian distributions"""
    def __init__(self, A, mu, sigma, covar_matrix):
        self.A = A
        self.mu = mu
        self.sigma = sigma
        self.covar_matrix = covar_matrix
    def value(self, x):
        return gaussfcn(x, self.A, self.mu, self.sigma)
    def area(self):
        return np.sqrt(2*np.pi)*self.A*np.abs(self.sigma)
    def as_string(self, ndigits=4):
        return str("A: {}, mu: {}, sigma: {}".format(round(self.A, ndigits),
                                                     round(self.mu, ndigits),
                                                     round(self.sigma, ndigits)))
    def print_full_info(self):
        print("Estimates of (A mu sigma) = (", self.A, self.mu, self.sigma, ")\n")
        print("Covariance matrix = \n", self.covar_matrix, "\n")
        print("Uncertainties in the estimated parameters: \n[ sigma^2(A) sigma^2(mu), sigma^2(sigma) ] = \n[", self.covar_matrix[0][0], self.covar_matrix[1][1], self.covar_matrix[2][2], "]\n" )
    

def fit_Gaussian(x, y, mu_guess, n):
    """ a simple function that tries to fit a Gaussian and return a Gauss object if fit was successful """
    peak_center_idx = (np.abs(x-mu_guess)).argmin() # find index of the mu guess value in bin_centers array 
    A_guess = y[peak_center_idx]                                   # a guess for the amplitude of the peak 
    sigma_guess = 1                                                # guess for sigma 
    guess = [A_guess, mu_guess, sigma_guess]
    
    # select values from bin_centers and y arrays of your data that correspond to your initiall guess af a peak
    peak_x = x[peak_center_idx-n:peak_center_idx+n]
    peak_y = y[peak_center_idx-n:peak_center_idx+n]

    ## scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
    import warnings
    from scipy.optimize import OptimizeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        ## perform the gaussian fit to the data:
        try:
            ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
            ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
            estimates, covar_matrix = curve_fit(GaussFunc,
                                    peak_x,
                                    peak_y,
                                    p0=guess)
            ## create a Gauss object with the fitted coefficients for better code readability
            g_final = Gauss(estimates[0], estimates[1], estimates[2], covar_matrix)
            return g_final
        except (RuntimeError, OptimizeWarning, TypeError):
            print("Gaussian fit failed! Try specifying another mu_guess.")
            return 0


def perform_Gaussian_fit(x,y, mu_guess, n, left_selection=None, right_selection=None, plotting = 1, printing = 1):
    
    peak_center_idx = (np.abs(x-mu_guess)).argmin() # find index of the mu guess value in bin_centers array 

    # select values from bin_centers and y arrays of your data that correspond to your initiall guess af a peak
    peak_x = x[peak_center_idx-n:peak_center_idx+n]
    peak_y = y[peak_center_idx-n:peak_center_idx+n]
    
  
    if (left_selection and right_selection):
        ############ Selecting points to fit linear function 
        #find indices of values in x array
        left_idx = [(np.abs(x-left_selection[0])).argmin(), (np.abs(x-left_selection[1])).argmin()]
        right_idx = [(np.abs(x-right_selection[0])).argmin(), (np.abs(x-right_selection[1])).argmin()]
        
        # select parts of data 
        left_x = x[left_idx[0]:(left_idx[1]+1)]
        right_x = x[right_idx[0]:(right_idx[1]+1)]
        left_y = y[left_idx[0]:(left_idx[1]+1)]
        right_y = y[right_idx[0]:(right_idx[1]+1)]
        
        # join left and right data points into one array for linear fitting 
        lin_x = np.concatenate([left_x, right_x])
        lin_y = np.concatenate([left_y, right_y])

        ############ Fitting linear function to selected points 
        guess = [2, 1] # guess parameters for linear fit 
        estimates_lin, covar_matrix = curve_fit(LineFunc,
                                            lin_x,
                                            lin_y,
                                            p0 = guess)

        #print("Linear fit coefficients (k m) = (", estimates_lin[0], estimates_lin[1], ")\n")
        ############ Subtract area under the linear fit
        peak_lin = LineFunc(peak_x, estimates_lin[0], estimates_lin[1])
        y_subst = peak_y - peak_lin
        peak_y = y_subst.copy()    
        
    ############ Fit a Gaussian to the peak without background    
    g_final = fit_Gaussian(peak_x, peak_y, mu_guess, n)   
    
    # Check if fit worked 
    if(g_final==0):
           return
       
    if (plotting):
        plt.figure()
        
        #Choose colour of the Gaussian depending on if fit was okay
        if (g_final.covar_matrix[2][2]<g_final.sigma):
            color = 'forestgreen'
        else:
            color = 'r'
        
        #plotting data
        plt.step(peak_x, y[peak_center_idx-n:peak_center_idx+n], where='mid', color='cornflowerblue', label='data')
        
        if (left_selection and right_selection):
            #plot points to which linear function is fitted 
            plt.step(left_x, left_y, where='mid', color='y')
            plt.step(right_x, right_y, where='mid', color='y')
            #plot support lines around selected points 
            plt.plot([left_idx[0]+0.5, left_idx[0]+0.5], [y[left_idx[0]]+0.5, g_final.A], color='y', linestyle="--")
            plt.plot([left_idx[1]+0.5, left_idx[1]+0.5], [y[left_idx[1]]+0.5, g_final.A], color='y', linestyle="--")
            plt.plot([right_idx[0]+0.5, right_idx[0]+0.5], [y[right_idx[0]]+0.5, g_final.A], color='y', linestyle="--")
            plt.plot([right_idx[1]+0.5, right_idx[1]+0.5], [y[right_idx[1]]+0.5, g_final.A], color='y', linestyle="--")
            # plot Gaussian 
            plt.plot(peak_x, peak_lin + GaussFunc(peak_x, g_final.A, g_final.mu, g_final.sigma), color=color, label = 'Gaussian fit')
            # plot linear fit
            plt.plot(lin_x, LineFunc(lin_x, estimates_lin[0], estimates_lin[1]), color='r', label = 'linear fit', alpha=0.6) 
        else:
            plt.plot(peak_x, GaussFunc(peak_x, g_final.A, g_final.mu, g_final.sigma), color=color, label = 'Gaussian fit')
        plt.legend(loc='upper right', frameon=False)
        plt.show()
        
    if (printing):
        g_final.print_full_info()
        
    return g_final
