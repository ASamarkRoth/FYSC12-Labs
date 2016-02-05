import csv                               # for reading in our data files
import logging                           # for orderly print output
import numpy as np                       # for handling of data
import matplotlib.pyplot as plt          # for plotting
from scipy.optimize import curve_fit     # to fit functions to the data
from scipy import signal                 # to find maxima in the data

class Spectrum:
    """A class to hold our spectrum measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(1))    # creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(1))    # creates a new empty array; will later store our y values
        self.name = filename   # a more descriptive name, can be used e.g. in legends
        self.duration = 0
    def subtract_from_data(self, m):
        self.y = self.y - m.y
        # these are spectra: cannot have counts below 0, so remove these here and set them to 0 instead:
        self.y[self.y < 0] = 0;
    def scale_data(self, scale):
        self.y *= scale

def gauss(x, *p):
    """ gauss function to be used for fits to the data"""
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

        
def read_mca_data_file(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Spectrum' object. Tested with Maestro Version 6.05 """
    log = logging.getLogger('betalab_analysis') # set up logging
    m = Spectrum(filename) # create a new Spectrum measurement object; this is what we return in the end
    log.info("Reading data from file '" + filename + "'")
    with open(filename, newline='') as f:
        reader = csv.reader(f) # use the python csv module to parse the file
        interval = []          # start/stop channel numbers used to assign correct x values to the data points
        # first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
        for row in reader:
            if row[0] == '$MEAS_TIM:':
                # this item gives the duration of the measurement
                log.debug("Parsing MEAS_TIM header info")
                row = next(reader)
                duration = [int(s) for s in row[0].split(' ')]
                m.duration = duration[1] # two parts: CPU/realtime; take the second
            if row[0] == '$DATA:':
                # this is the last part of the header and contains the start/stop channel numbers
                log.debug("Parsing DATA header info")
                row = next(reader)
                interval = [int(s) for s in row[0].split(' ')]
                # "DATA" is the last item: stop with the header processing
                break
        # TODO!!! make sure that the file does not end before we have parsed the header!
        log.debug("Done with header parsing")
        nbins = int(interval[1]-interval[0])+1
        m.y = np.array(np.zeros(nbins))
        # continue, now reading data
        for idx, row in enumerate(reader):
            if idx >= nbins:
                break
            m.y[idx] = int(row[0])
        m.x = np.arange(interval[0], interval[1]+1,1)
        log.debug("Loaded all data from file")
    return m

def fit_gaussians_to_measurement(m):
    """ fits all gaussians in a spectrum measurement and returns a list of coefficients"""
    # list to store the paramters of the fitted gaussians in
    gaussians = []
    
    # find peaks in m.y with widths given by an array
    peakind = signal.find_peaks_cwt(m.y, np.arange(10,80,5)) 

    for p in peakind:
        log.info("Found peak in the data at position x = " + str(m.x[p]) + "; fitting it with a gaussian")
        # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
        p0 = [1., m.x[p], 1.] # mu is given by one of the found peaks positions
        # scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
        import warnings
        from scipy.optimize import OptimizeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            # perform the gaussian fit to the data:
            try:
                #coeff, var_matrix = curve_fit(gauss, m.x, m.y, p0=p0) # fit using the full data range
                coeff, var_matrix = curve_fit(gauss, m.x[p-10:p+10], m.y[p-10:p+10], p0=p0) # fit using "slices" of the arrays with +/- 10 around peak
            except (RuntimeError, OptimizeWarning):
                # the minimization did not work out... log it and continue to next peak
                log.info("  - gaussian fit failed!")
                continue
        # filter the results
        xdynrange = m.x.shape[0] # dynamic range in x
        if coeff[2] > 0.1*xdynrange: # check width of gaussian in percent of the dynamic range
            log.info("  - sigma out of bounds: " + str(coeff[2]))
            continue
        if coeff[1] > xdynrange-0.1*xdynrange or coeff[1] < 6: # check center position: should not be at the limits of measurement range
            log.info("  - mu out of bounds: " + str(coeff[1]))
            continue
        if coeff[0] < 0.25*np.average(m.y): # check if the peak hight is at least 25% over the average data value
            log.info("  - A out of bounds: " + str(coeff[0]))
            continue            
        log.info("  - fit result: A = " + str(coeff[0]) + ", mu = " + str(coeff[1]) + ", sigma = " + str(coeff[2]) + ". ")
        # store the results
        gaussians.append(coeff)
    return gaussians
        
if __name__ == '__main__':
    # set up some print-out routines (logging)
    FORMAT = '%(asctime)s %(name)s %(levelname)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('betalab_analysis') # set up logging
    log.setLevel("DEBUG")

    # read in the first file
    m = read_mca_data_file('data samples/Beta NY 2015/P32 60 min.Spe')

    # make the plot pretty
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title(m.name)
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$') # to add text to the plot
    # plt.axis([40, 160, 0, 0.03]) # to set the axis range
    # plt.yscale('log') # set y axis to log scale
    plt.grid(True)

    # plot the first file
    plt.plot(m.x, m.y, 'o')       # plot with markers

    # read in the Cs137 file measured without Al plate
    cs137 = read_mca_data_file('data samples/Beta NY 2015/cs137 15 min utan Al.Spe')

    # plot into a new figure
    plt.figure(2)
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title("Cs-137")

    plt.plot(cs137.x, cs137.y, 'o',label="Cs-137")       # plot with markers

    # read in the gamma background measurement of Cs137
    cs137_gamma = read_mca_data_file('data samples/Beta NY 2015/cs137 15 min med Al.Spe')

    # weight the spectrum according to the ratio of measurement durations
    cs137_gamma.scale_data(cs137.duration/cs137_gamma.duration)
    
    # plot into same figure
    plt.plot(cs137_gamma.x, cs137_gamma.y, 'o', label="Cs-137 gamma bkgrd")

    # now subtract the background from the measurement
    cs137.subtract_from_data(cs137_gamma)

    # plot the resut
    plt.plot(cs137.x, cs137.y, 'o', label="Cs-137 w/o gamma bkgrd")

    # fit all gaussians in our measurement
    fits = fit_gaussians_to_measurement(cs137)

    # loop over fit results
    for g in fits:
        # plot the gaussian fit
        plt.plot(m.x, gauss(cs137.x,*g), label="Gauss fit, $\sigma$="+str(g[2]))

    
    # generate the legend (with the "label" information from the plots)
    plt.legend()

    # now we have some data for our energy calibration:
    # peak 0: pedestal, energy 0
    # peak 1: internal convertion peak, 0.630 MeV
    ecalib_data_cs137 = np.array([fits[0][1], 0 ],
                                 fits[1][1] , 0.630 )

    # plot into a new figure
    plt.figure(3)
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('energy [MeV]')
    plt.title("Energy Calibration")

    plt.plot(ecalib_data_cs137, 'o',label="Cs-137")

    
    # final step:
    plt.show()           # <-- shows the plot (not needed with interactive plots)
