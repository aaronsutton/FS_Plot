import numpy as np
from scipy.interpolate import griddata


def return_data(filename):
    f_plot = open('data_' + filename + '.dat', 'r')
    plotData = np.loadtxt(f_plot)
    plotX = plotData[:, 0]
    plotY = plotData[:, 1]
    plotZ = plotData[:, 2]
    Energy = plotData[:, 3]
    R = plotData[:, 0:3]
    pts = 30j
    X, Y, Z = np.mgrid[min(plotX):max(plotX):pts, min(plotY):max(plotY):pts, min(plotZ):max(plotZ):pts]
    F = griddata(R, Energy, (X, Y, Z), method='linear')
    f_plot.close()
    return X, Y, Z, F
