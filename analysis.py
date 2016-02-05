import csv
import logging
import numpy as np
import matplotlib.pyplot as plt  # the tidy way

class Measurement:
    """A class to hold our measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(1))    # creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(1))    # creates a new empty array; will later store our y values
        self.name = filename   # a more descriptive name, can be used e.g. in legends
        self.duration = 0

def readMcaDataFile(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Measurement' object. Tested with Maestro Version 6.05 """
    log = logging.getLogger('betalab_analysis') # set up logging
    m = Measurement(filename) # create a new Measurement object; this is what we return in the end
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
        log.debug("Loaded data from file")
        return m
        
if __name__ == '__main__':
    # set up some print-out routines (logging)
    log = logging.getLogger('betalab_analysis') # set up logging
    formatter = logging.Formatter('%(asctime)s %(name)s(%(levelname)s): %(message)s',"%H:%M:%S")
    log.setLevel("DEBUG")

    # read in the first file
    m = readMcaDataFile('data samples/Beta NY 2015/P32 60 min.Spe')
    
    # plot the first file
    plt.plot(m.x, m.y, 'o')       # plot with markers

    # make the plot pretty
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title(m.name)
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$') # to add text to the plot
    # plt.axis([40, 160, 0, 0.03]) # to set the axis range
    plt.yscale('log')
    plt.grid(True)

    # final step:
    plt.show()           # <-- shows the plot (not needed with interactive plots)
