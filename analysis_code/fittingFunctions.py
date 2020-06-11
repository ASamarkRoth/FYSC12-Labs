#!/usr/bin/env python3
import csv                               ## for reading in our data files
import logging                           ## for orderly print output
import numpy as np                       ## for handling of data
import sys                               ## useful system calls (used to exit cleanly)
import os                                ## path manipulations
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
    

def perform_Gaussian_fit(data_x, data_y, mu_guess, n, left_selection=None, right_selection=None, plotting_main = 1, printing = 1, plotting_subt=1):
    # make a copy of original data:
    x = data_x.copy()
    y = data_y.copy()
    
    """Fit a Guassian to data"""
    if (left_selection and right_selection):
        
        ############ Selecting points to fit linear function 
        x_selected = np.concatenate([x[left_selection[0]:(left_selection[1]+1)], x[right_selection[0]:(right_selection[1]+1)]])
        y_selected = np.concatenate([y[left_selection[0]:(left_selection[1]+1)], y[right_selection[0]:(right_selection[1]+1)]])

        if (plotting_subt): 
            print("Selected data regions to fit the line:")
            plt.figure()  
            #plotting data
            plt.step(data_x[mu_guess-n:mu_guess+n], data_y[mu_guess-n:mu_guess+n], where='mid')
            #plot points to which linear function is fitted 
            plt.step(x[left_selection[0]:(left_selection[1]+1)], y[left_selection[0]:(left_selection[1]+1)], where='mid', color='y')
            plt.step(x[right_selection[0]:(right_selection[1]+1)], y[right_selection[0]:(right_selection[1]+1)], where='mid', color='y')
            #plot dashed lines  
            plt.plot([left_selection[0]-0.5, left_selection[0]-0.5], [y[left_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([left_selection[1]+0.5, left_selection[1]+0.5], [y[left_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[0]-0.5, right_selection[0]-0.5], [y[right_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[1]+0.5, right_selection[1]+0.5], [y[right_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.show()

        ############ Fitting linear function to selected points 
        guess = [2, 1]
        estimates, covar_matrix = curve_fit(LineFunc,
                                            x_selected,
                                            y_selected,
                                            p0 = guess)

        print("Linear fit coefficients (k m) = (", estimates[0], estimates[1], ")\n")
        
        if(plotting_subt):
            plt.figure()
            #plotting data
            plt.step(data_x[mu_guess-n:mu_guess+n], data_y[mu_guess-n:mu_guess+n], where='mid')
            #plot points to which linear function is fitted 
            plt.step(x[left_selection[0]:(left_selection[1]+1)], y[left_selection[0]:(left_selection[1]+1)], where='mid', color='y')
            plt.step(x[right_selection[0]:(right_selection[1]+1)], y[right_selection[0]:(right_selection[1]+1)], where='mid', color='y')
            #plot dashed lines  
            plt.plot([left_selection[0]-0.5, left_selection[0]-0.5], [y[left_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([left_selection[1]+0.5, left_selection[1]+0.5], [y[left_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[0]-0.5, right_selection[0]-0.5], [y[right_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[1]+0.5, right_selection[1]+0.5], [y[right_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            # plot linear fit
            plt.plot(x_selected, LineFunc(x_selected, estimates[0], estimates[1]), color='r', label='linear fit') 
            plt.show()

        ############ Subtract area under the linear fit

        y_lin = LineFunc(x[mu_guess-n:mu_guess+n], estimates[0], estimates[1])
        y_before =  y[mu_guess-n:mu_guess+n].copy() # just for plotting 
        y[mu_guess-n:mu_guess+n] = y[mu_guess-n:mu_guess+n] - y_lin

        if (plotting_subt):
            plt.figure()
            #plotting data
            plt.step(x[mu_guess-n:mu_guess+n], data_y[mu_guess-n:mu_guess+n], where='mid', label = 'original data')
            #plot points to which linear function is fitted 
            plt.step(x[left_selection[0]:(left_selection[1]+1)], y[left_selection[0]:(left_selection[1]+1)], where='mid', color='y')
            plt.step(x[right_selection[0]:(right_selection[1]+1)], y[right_selection[0]:(right_selection[1]+1)], where='mid', color='y')
            #plot dashed lines  
            plt.plot([left_selection[0]-0.5, left_selection[0]-0.5], [y[left_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([left_selection[1]+0.5, left_selection[1]+0.5], [y[left_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[0]-0.5, right_selection[0]-0.5], [y[right_selection[0]]+0.5, y[mu_guess]], color='y', linestyle="--")
            plt.plot([right_selection[1]+0.5, right_selection[1]+0.5], [y[right_selection[1]]+0.5, y[mu_guess]], color='y', linestyle="--")
            # plot linear fit
            plt.plot(x_selected, LineFunc(x_selected, estimates[0], estimates[1]), color='r', label='linear fit') 
            # plot data with substracted background 
            plt.step(x[mu_guess-n:mu_guess+n], y[mu_guess-n:mu_guess+n], where='mid', color='g')
            plt.legend(loc='best', frameon=False)
            plt.show()

        
    ############ Fit a Gaussian to data ######################################   

    A_guess = y[mu_guess]
    sigma_guess = 1
    guess = [A_guess, mu_guess, sigma_guess]

    estimates, covar_matrix = curve_fit(GaussFunc,
                                    x,
                                    y,
                                    p0=guess)
    g_final = Gauss(estimates[0], estimates[1], estimates[2], covar_matrix )
    
    color = "forestgreen"
#     if ((1-mu_guess/g_final.mu) <= 0.005):
#         color = 'forestgreen'
#     else:
#         color = 'r'
           
    
    if (plotting_main):
        plt.figure()
        plt.step(x[mu_guess-n:mu_guess+n], y[mu_guess-n:mu_guess+n], where='mid', color='cornflowerblue', label='data')
        plt.plot(x[mu_guess-n:mu_guess+n], GaussFunc(x[mu_guess-n:mu_guess+n], g_final.A, g_final.mu, g_final.sigma), color=color, label = 'Gaussian fit')
        plt.legend(loc='upper right', frameon=False)
        plt.show()
        
    if (printing):
        print("Estimates of (A mu sigma) = (", g_final.A, g_final.mu, g_final.sigma, ")\n")
        print("Covariance matrix = \n", g_final.covar_matrix, "\n")
        print("Uncertainties in the estimated parameters: \n[ sigma^2(A) sigma^2(mu), sigma^2(sigma) ] = \n[", g_final.covar_matrix[0][0], g_final.covar_matrix[1][1], g_final.covar_matrix[2][2], "]\n" )
 
    return g_final