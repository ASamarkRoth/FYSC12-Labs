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
    def as_string(self, ndigits=4):
        return str("A: {}, mu: {}, sigma: {}".format(round(self.A, ndigits),
                                                     round(self.mu, ndigits),
                                                     round(self.sigma, ndigits)))
    

def fit_Gaussian(x, y, mu_guess, n):
    """ a simple function that tries to fit a Gaussian and return a Gauss object if fit was successful """
    idx = (np.abs(x-mu_guess)).argmin()
    A_guess = y[idx]
    sigma_guess = 1
    guess = [A_guess, mu_guess, sigma_guess]
    
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
                                    x[idx-n:idx+n],
                                    y[idx-n:idx+n],
                                    p0=guess)
            ## create a Gauss object with the fitted coefficients for better code readability
            g_final = Gauss(estimates[0], estimates[1], estimates[2], covar_matrix)
            return g_final
        except (RuntimeError, OptimizeWarning, TypeError):
            print("Gaussian fit failed! Try specifying another mu_guess.")
            return 0


def perform_Gaussian_fit(data_x, data_y, mu_guess, n, left_selection=None, right_selection=None, plotting_main = 1, printing = 1):
    x = data_x.copy()
    y = data_y.copy()
    idx = (np.abs(x-mu_guess)).argmin()
    A_guess = y[idx]
    sigma_guess = 1
    guess = [A_guess, mu_guess, sigma_guess]
    
    if (left_selection and right_selection):
        ############ Selecting points to fit linear function 
        x_selected = np.concatenate([x[left_selection[0]:(left_selection[1]+1)], x[right_selection[0]:(right_selection[1]+1)]])
        y_selected = np.concatenate([y[left_selection[0]:(left_selection[1]+1)], y[right_selection[0]:(right_selection[1]+1)]])

        ############ Fitting linear function to selected points 
        guess = [2, 1]
        estimates_lin, covar_matrix = curve_fit(LineFunc,
                                            x_selected,
                                            y_selected,
                                            p0 = guess)

        print("Linear fit coefficients (k m) = (", estimates_lin[0], estimates_lin[1], ")\n")
        ############ Subtract area under the linear fit
        y_lin = LineFunc(x[idx-n:idx+n], estimates_lin[0], estimates_lin[1])
        y_subst = y[idx-n:idx+n] - y_lin
        y[idx-n:idx+n] = y_subst
            
    ############ Fit a Gaussian to data ######################################   
    g_final = fit_Gaussian(x, y, mu_guess, n)
    
    if(g_final==0):
           return

    if (g_final.covar_matrix[2][2]<g_final.sigma):
        color = 'forestgreen'
    else:
        color = 'r'
       
    if (plotting_main):
        plt.figure()
        #plotting data
        plt.step(data_x[idx-n:idx+n], data_y[idx-n:idx+n], where='mid', color='cornflowerblue', label='data')
        
        if (left_selection and right_selection):
            #plot points to which linear function is fitted 
            plt.step(data_x[left_selection[0]:(left_selection[1]+1)], data_y[left_selection[0]:(left_selection[1]+1)], where='mid', color='y')
            plt.step(data_x[right_selection[0]:(right_selection[1]+1)], data_y[right_selection[0]:(right_selection[1]+1)], where='mid', color='y')
            #plot support lines  
            plt.plot([left_selection[0]-0.5, left_selection[0]-0.5], [data_y[left_selection[0]]+0.5, data_y[idx]], color='y', linestyle="--")
            plt.plot([left_selection[1]+0.5, left_selection[1]+0.5], [data_y[left_selection[1]]+0.5, data_y[idx]], color='y', linestyle="--")
            plt.plot([right_selection[0]-0.5, right_selection[0]-0.5], [data_y[right_selection[0]]+0.5, data_y[idx]], color='y', linestyle="--")
            plt.plot([right_selection[1]+0.5, right_selection[1]+0.5], [data_y[right_selection[1]]+0.5, data_y[idx]], color='y', linestyle="--")
            # plot Gaussian 
            plt.plot(x[idx-n:idx+n], y_lin + GaussFunc(x[idx-n:idx+n], g_final.A, g_final.mu, g_final.sigma), color=color, label = 'Gaussian fit')
            # plot linear fit
            plt.plot(x_selected, LineFunc(x_selected, estimates_lin[0], estimates_lin[1]), color='r', label = 'linear fit', alpha=0.6) 
        else:
            plt.plot(x[idx-n:idx+n], GaussFunc(x[idx-n:idx+n], g_final.A, g_final.mu, g_final.sigma), color=color, label = 'Gaussian fit')
        plt.legend(loc='upper right', frameon=False)
        plt.show()
        
    if (printing):
        print("Estimates of (A mu sigma) = (", g_final.A, g_final.mu, g_final.sigma, ")\n")
        print("Covariance matrix = \n", g_final.covar_matrix, "\n")
        print("Uncertainties in the estimated parameters: \n[ sigma^2(A) sigma^2(mu), sigma^2(sigma) ] = \n[", g_final.covar_matrix[0][0], g_final.covar_matrix[1][1], g_final.covar_matrix[2][2], "]\n" )
 
    return g_final

# def plot_interactive_fit(data_xx, y, mu_guess, n):
#     gauss = perform_Gaussian_fit(x, y, mu_guess, n)

# def perform_Gaussian_fit_with_widget(x, y, mu_guess, n):
#     interactive_plot = interactive(plot_interactive_fit, x=fixed(x), data_y=fixed(y), mu_guess=(3200, 3400, 1), n=(15, 45, 1), continuous_update=False)
#     interactive_plot.children[0].description=r'mu_guess' # slider
#     interactive_plot.children[1].description=r'n'
#     interactive_plot.children[0].continuous_update = True
#     interactive_plot.children[1].continuous_update = False
#     return interactive_plot
