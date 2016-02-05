import csv
import logging
import numpy as np
import matplotlib.pyplot as plt  # the tidy way

class Measurement:
    """A class to hold our measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array()    # creates a new empty array; will later store our x values
        self.y = np.array()    # creates a new empty array; will later store our y values
        self.name = filename   # a more descriptive name, can be used e.g. in legends
        self.duration = 0

def readMcaDataFile(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Measurement' object. Tested with Maestro Version 6.05 """
    interval = []
    datay = np.array(np.zeros(512))
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        # first parse the header of the data file (until the '$DATA:' line)
        for row in reader:
            if row[0] == '$DATA:':
                row = next(reader)
                interval = [int(s) for s in row[0].split(' ')]
                print("found interval")
                # stop with the header processing
                break
        # !!! make sure that the file does not end before we have parsed the header!
        nbins = int(interval[1]-interval[0])
        np.resize(datay,nbins)
        # continue, now reading data
        for idx, row in enumerate(reader):
            if idx > nbins:
                break
            datay[idx] = int(row[0])
    x = np.arange(interval[0], interval[1]+1,1)
	m = Measurement(filename) # create a new Measurement object

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)  
    readMcaDataFile('data samples/Beta NY 2015/P32 60 min.Spe')
    plt.plot(x, datay, fmt='o')       # plot with markers
    plt.show()           # <-- shows the plot (not needed with interactive plots)
