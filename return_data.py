import numpy as np
from scipy.interpolate import griddata


def return_data(filename, G_1, G_2, G_3):
    f_plot = open('data_' + filename + '.dat', 'r')
    plotData = np.loadtxt(f_plot)
    plotX = plotData[:, 0]
    plotY = plotData[:, 1]
    plotZ = plotData[:, 2]
    plotXorg = plotX.copy()
    plotYorg = plotY.copy()
    plotZorg = plotZ.copy()
    Energy = plotData[:, 3]
    Energyorg = Energy.copy()
    R = plotData[:, 0:3]

    shift_it = 1
    shift_again = 0
    if shift_it:
        translate = [1, 1, -1]
        shift = translate[0] * np.array(G_1) + translate[1] * np.array(G_2) + translate[2] * np.array(G_3)
        shift_x = [x + shift[0] for x in plotX]
        shift_y = [x + shift[1] for y in plotY]
        shift_z = [x + shift[2] for z in plotZ]
        plotX = np.append(plotX, shift_x)
        plotY = np.append(plotY, shift_y)
        plotZ = np.append(plotZ, shift_z)
        newR = zip(shift_x, shift_y, shift_z)
        R = np.append(R, newR, axis=0)
        Energy = np.append(Energy, Energyorg)

        if shift_again:
            Energy = np.append(Energy, Energyorg)
            translate = [0, 1, 0]
            shift = translate[0] * G_1 + translate[1] * G_2 + translate[2] * G_3
            shift_x = [x + shift[0] for x in plotXorg]
            shift_y = [x + shift[1] for y in plotYorg]
            shift_z = [x + shift[2] for z in plotZorg]
            plotX = np.append(plotX, shift_x)
            plotY = np.append(plotY, shift_y)
            plotZ = np.append(plotZ, shift_z)
            newR2 = zip(shift_x, shift_y, shift_z)
            R = np.append(R, newR2, axis=0)

    pts = 100j
    X, Y, Z = np.mgrid[min(plotX):max(plotX):pts, min(plotY):max(plotY):pts, min(plotZ):max(plotZ):pts]
    # X, Y, Z = np.mgrid[min(bzx):max(bzx):pts, min(bzy):max(bzy):pts, min(bzz):max(bzz):pts]
    F = griddata(R, Energy, (X, Y, Z), method='linear')
    f_plot.close()
    return X, Y, Z, R, F
