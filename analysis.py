#!/usr/bin/env python3
import csv                               ## for reading in our data files
import logging                           ## for orderly print output
import numpy as np                       ## for handling of data
import matplotlib.pyplot as plt          ## for plotting
from scipy.optimize import curve_fit     ## to fit functions to the data
from scipy import signal                 ## to find maxima in the data
from scipy import stats                  ## for linear regressions
import sys                               ## useful system calls (used to exit cleanly)

class Spectrum:
    """A class to hold our spectrum measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(1))    ## creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(1))    ## creates a new empty array; will later store our y values
        self.name = filename   ## a more descriptive name, can be used e.g. in legends
        self.duration = 0
    def subtract_from_data(self, m):
        self.y = self.y - m.y
        ## these are spectra: cannot have counts below 0, so remove these here and set them to 0 instead:
        self.y[self.y < 0] = 0;
    def scale_data(self, scale):
        self.y *= scale
    def energy_calibrate(self, slope, intercept):
        self.x = self.x*slope + intercept

def gauss(x, *p):
    """ gauss function to be used for fits to the data"""
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def line(x, *p):
    """ straight line function to be used for fits to the data"""
    a, b = p
    return a*x+b

        
def read_mca_data_file(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Spectrum' object. Tested with Maestro Version 6.05 """
    log = logging.getLogger('betalab_analysis') ## set up logging
    m = Spectrum(filename) ## create a new Spectrum measurement object; this is what we return in the end
    log.info("Reading data from file '" + filename + "'")
    try:
        with open(filename, newline='') as f:
            reader = csv.reader(f) ## use the python csv module to parse the file
            interval = []          ## start/stop channel numbers used to assign correct x values to the data points
            ## first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
            for row in reader:
                if row[0] == '$MEAS_TIM:':
                    ## this item gives the duration of the measurement
                    log.debug("Parsing MEAS_TIM header info")
                    row = next(reader)
                    duration = [int(s) for s in row[0].split(' ')]
                    m.duration = duration[1] ## two parts: CPU/realtime; take the second
                if row[0] == '$DATA:':
                    ## this is the last part of the header and contains the start/stop channel numbers
                    log.debug("Parsing DATA header info")
                    row = next(reader)
                    interval = [int(s) for s in row[0].split(' ')]
                    ## "DATA" is the last item: stop with the header processing
                    break
            ## TODO!!! make sure that the file does not end before we have parsed the header!
            log.debug("Done with header parsing")
            nbins = int(interval[1]-interval[0])+1
            m.y = np.array(np.zeros(nbins))
            ## continue, now reading data
            for idx, row in enumerate(reader):
                if idx >= nbins:
                    break
                m.y[idx] = int(row[0])
            m.x = np.arange(interval[0], interval[1]+1,1)
            log.debug("Loaded all data from file")
    except IOError:
        log.error("Could not find the file '"+str(filename)+"'")
        sys.exit(-1)
    return m

def fit_gaussians_to_measurement(m):
    """ fits all gaussians in a spectrum measurement and returns a list of coefficients"""
    ## list to store the paramters of the fitted gaussians in
    gaussians = []
    
    ## find peaks in m.y with widths given by an array
    peakind = signal.find_peaks_cwt(m.y, np.arange(10,80,5)) 

    for p in peakind:
        log.info("Found peak in the data at position x = " + str(m.x[p]) + "; fitting it with a gaussian")
        ## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
        p0 = [1., m.x[p], 1.] ## mu is given by one of the found peaks positions
        ## scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
        import warnings
        from scipy.optimize import OptimizeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            ## perform the gaussian fit to the data:
            try:
                ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
                ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
                ## coeff, var_matrix = curve_fit(gauss, m.x, m.y, p0=p0) ## fit using the full data range, might not work with multiple peaks
                coeff, var_matrix = curve_fit(gauss, m.x[p-10:p+10], m.y[p-10:p+10], p0=p0) ## fit using "slices" of the arrays with +/- 10 around peak
            except (RuntimeError, OptimizeWarning):
                ## the minimization did not work out... log it and continue to next peak
                log.info("  - gaussian fit failed!")
                continue
        ## filter the results -- or we get a lot of "spurious" peaks
        xdynrange = m.x.shape[0] ## dynamic range in x
        if coeff[2] > 0.1*xdynrange: ## check width of gaussian in percent of the dynamic range
            log.info("  - sigma out of bounds: " + str(coeff[2]))
            continue
        if coeff[1] > xdynrange-0.1*xdynrange or coeff[1] < 6: ## check center position: should not be at the limits of measurement range
            log.info("  - mu out of bounds: " + str(coeff[1]))
            continue
        if coeff[0] < 0.25*np.average(m.y): ## check if the peak hight is at least 25% over the average data value
            log.info("  - A out of bounds: " + str(coeff[0]))
            continue            
        log.info("  - fit result: A = " + str(coeff[0]) + ", mu = " + str(coeff[1]) + ", sigma = " + str(coeff[2]) + ". ")
        ## store the results
        gaussians.append(coeff)
    return gaussians
        
if __name__ == '__main__':
    ## set up some print-out routines (logging)
    FORMAT = '%(asctime)s %(name)s %(levelname)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('betalab_analysis') ## set up logging
    log.setLevel("INFO")

    ##            _________
    ## _ __      |___ /___ \
    ##| '_ \ _____ |_ \ __) |
    ##| |_) |_____|__) / __/
    ##| .__/     |____/_____|
    ##|_|    
    ## setup the plot
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title("P-32 MCA spectrum")
    ## plt.text(60, .025, r'$\mu=100,\ \sigma=15$') ## to add text to the plot
    ## plt.axis([40, 160, 0, 0.03]) ## to set the axis range
    plt.yscale('log') ## set y axis to log scale
    plt.grid(True)

    ## read in the p32 measurement file
    p32 = read_mca_data_file('data/p32.Spe')
    ## plot the p32 raw measurement
    plt.plot(p32.x, p32.y, 'o', label="P-32 raw data")  ## 'o' parameter: plot with markers

    
    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...

    
    ## read in the background measurement file
    bkgrd = read_mca_data_file('data/p32-background.Spe')
    ## normalize the background to the measurement time ratio between P32 and bkgrd
    log.info("Adjusting P-32 background normalization by: " + str(p32.duration/bkgrd.duration))
    bkgrd.scale_data(p32.duration/bkgrd.duration)
    ## show background also in plot
    plt.plot(bkgrd.x, bkgrd.y, 'o', label="background")

    ## now subtract the background from the p32 measurement
    p32.subtract_from_data(bkgrd)
    ## plot the background-subracted p32 measurement
    plt.plot(p32.x, p32.y, 'o', label="P-32 - bkgrd")  ## 'o' parameter: plot with markers
    plt.legend()     ## generate the legend (with the "label" information from the plots)

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
        
    
    ##                _ __________
    ##  ___ ___      / |___ /___  |
    ## / __/ __|_____| | |_ \  / /
    ##| (__\__ \_____| |___) |/ /
    ## \___|___/     |_|____//_/
    ##    
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title("Cs-137")

    ## read in the Cs137 file measured without Al plate
    cs137 = read_mca_data_file('data/cs137.Spe')
    ## plot the cs-137 data
    plt.plot(cs137.x, cs137.y, 'o',label="Cs-137")

    ## read in the gamma background measurement of Cs137
    cs137_gamma = read_mca_data_file('data/cs137-background.Spe')
    ## weight the spectrum according to the ratio of measurement durations
    cs137_gamma.scale_data(cs137.duration/cs137_gamma.duration)
    ## plot the cs-137 e-suppressed data into same figure
    plt.plot(cs137_gamma.x, cs137_gamma.y, 'o', label="Cs-137 gamma bkgrd")

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    
    ## now subtract the background from the measurement
    cs137.subtract_from_data(cs137_gamma)
    ## plot the result
    plt.plot(cs137.x, cs137.y, 'o', label="Cs-137 w/o gamma bkgrd")

    ## fit all gaussians in our measurement
    fits = fit_gaussians_to_measurement(cs137)
    ## loop over fit results
    for g in fits:
        ## plot the gaussian fit
        plt.plot(cs137.x, gauss(cs137.x,*g), label="Gauss fit, $\sigma$="+str(g[2]))
    ## now we have some data for our energy calibration:
    ## peak 0: pedestal, energy 0
    ## peak 1: internal conversion peak, 0.630 MeV
    ecalib_data_cs137 = np.array([[fits[0][1], fits[1][1] ], #x
                                  [0. , 0.630 ]] ) ## y
    ## copy the sigma (the width of the gaussian) from the fit results
    ecalib_sigma_cs137 = np.array([fits[0][2], fits[1][2]])                                  
    ## generate the legend (with the "label" information from the plots)
    plt.legend()

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...

    
    ##                                               _ _ _               _   _                 
    ##  ___ _ __   ___ _ __ __ _ _   _      ___ __ _| (_) |__  _ __ __ _| |_(_) ___  _ __      
    ## / _ \ '_ \ / _ \ '__/ _` | | | |    / __/ _` | | | '_ \| '__/ _` | __| |/ _ \| '_ \     
    ##|  __/ | | |  __/ | | (_| | |_| |   | (_| (_| | | | |_) | | | (_| | |_| | (_) | | | |    
    ## \___|_| |_|\___|_|  \__, |\__, |    \___\__,_|_|_|_.__/|_|  \__,_|\__|_|\___/|_| |_|    
    ##                     |___/ |___/
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('energy [MeV]')
    plt.title("Energy Calibration")

    ## plot the data from Cs-137
    ## could use "plt.errorbar" to include uncertainties!
    plt.plot(ecalib_data_cs137[0], ecalib_data_cs137[1], 'o',label="Cs-137")

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    
    points_x = np.concatenate((ecalib_data_cs137[0])) ## concatenate: join array data together; useful when we have more data to actually join...
    points_y = np.concatenate((ecalib_data_cs137[1]))
    sigma    = np.concatenate((ecalib_sigma_cs137))

    ## linear regression of the data
    ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.linregress.html
    slope, intercept, r_value, p_value, std_err = stats.linregress(points_x, points_y)

    ## ALTERNATIVE METHODS TO FIT:
    ## p0 = [1., 1.] ## initial guess for the fitting (coefficients a, b)
    ## coeff, var_matrix = curve_fit(line, points_x, points_y, sigma=sigma, p0=p0)
    ## slope = coeff[0]
    ## intercept = coeff[1]

    log.info("Determined calibration constants from linear regression: E [MeV] = "+str(slope)+"*N_ch + " + str(intercept))
    x = np.arange(1,512)
    plt.plot(x,slope*x+intercept,label="linear regression")
    plt.legend()
    
    ## apply energy calibration
    log.info("Applying calibration constants")
    p32.energy_calibrate(slope,intercept)
    
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('energy [MeV]')
    plt.ylabel('counts')
    plt.title("P-32 energy spectrum")
    plt.plot(p32.x, p32.y, 'o')

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
  
    
    ## _____                   _       _  __          _        ____  _       _
    ##|  ___|__ _ __ _ __ ___ (_)     | |/ /   _ _ __(_) ___  |  _ \| | ___ | |_
    ##| |_ / _ \ '__| '_ ` _ \| |_____| ' / | | | '__| |/ _ \ | |_) | |/ _ \| __|
    ##|  _|  __/ |  | | | | | | |_____| . \ |_| | |  | |  __/ |  __/| | (_) | |_
    ##|_|  \___|_|  |_| |_| |_|_|     |_|\_\__,_|_|  |_|\___| |_|   |_|\___/ \__|    
    ## fermi-kurie calculations:
    mec2 = 0.510998910 ## MeV
    pc = np.sqrt((p32.x + mec2)**2 - mec2**2)
    A = pc/mec2
    f = 1.3604*A*A + 0.1973*A + 0.0439
    Ee = (p32.x + mec2)
    QminTe = np.sqrt((p32.y*pc)/(Ee*f))

    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('Te [MeV]')
    plt.ylabel('Q-Te')
    plt.title("P-32 Fermi-Kurie")
    plt.plot(p32.x, QminTe, 'o', label="data")

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
    log.info(" Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    
    ## linear regression of the FM plot
    ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.linregress.html
    ## the fit does not really work on the edges of the FM plot, so we take the region 0.2<E [MeV]<1.5
    lower_limit = np.where(p32.x>0.2)[0][0] ## first elements indicate first bin matching our criteria
    upper_limit = np.where(p32.x>1.5)[0][0]
    slope, intercept, r_value, p_value, std_err = stats.linregress(p32.x[lower_limit:upper_limit], QminTe[lower_limit:upper_limit])
    x = np.arange(0,2,0.05) ## generate x axis for fit result (start, stop, stepsize)
    plt.plot(x,slope*x+intercept,label="linear regression")
    plt.legend()

    ## now the Q value is determined by where the linear regression intersects with the x axis (Q-Te = 0)
    Q = -intercept/slope

    ## print results
    log.info("Determined linear regression to Fermi-Kurie plot: Q-Te = "+str(slope)+"*Te + " + str(intercept))
    log.info("===> Q value for P32: Q = "+str(Q)+" MeV ")
    plt.text(1.6, 25, 'Q = '+"{:.3f}".format(Q)+' MeV') ## to add text to the plot
    
    ## final step:
    plt.show()           ## <-- shows the plot (not needed with interactive plots)
