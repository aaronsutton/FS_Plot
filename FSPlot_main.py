#import re
import numpy as np
import matplotlib.pyplot as plt
import itertools
#import matplotlib.figure as fig
from mpl_toolkits.mplot3d import Axes3D

assert Axes3D

data = []
filename = 'YbSO25072013.bxsf.band-40'
with open('/home/asutton/Dropbox/WIEN2K/CASES/YbSO25072013/' + filename, 'r') as f_id:
    f_id.seek(0)
    line = f_id.readline()
    line = f_id.readline()
    E_F = float(line.split()[2])
    line = f_id.readline()
    line = f_id.readline()
    line = f_id.readline()
    line = f_id.readline()
    line = f_id.readline()
    no_bands = int(line.split()[0])
    print 'Number of bands is {}'.format(no_bands)
    line = f_id.readline()
    dimensions = int(line.split()[0]), int(line.split()[1]), int(line.split()[2])
    line = f_id.readline()
    line = f_id.readline()
    G_1 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
    line = f_id.readline()
    G_2 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
    line = f_id.readline()
    G_3 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
    line = f_id.readline()
    band = int(line.split()[1])
    print 'The band number is {}'.format(band)
    #print(G_1)
    # Now loop over and read the data into an array that is easy to format
    testdone = "END_BANDGRID_3D"
    flag = 0
    while not flag:
        line = f_id.readline()
        #print(line)
        if testdone in line:
            print('DONE!')
            flag = 1
        else:
            data += [float(value) for value in line.strip().split()]

data = np.array(data)
data2 = data.reshape([dimensions[0], dimensions[1], dimensions[2]])
K_x = 0
data2d = data2[K_x, :, :]
x = np.linspace(0, np.linalg.norm(G_3), dimensions[2])
y = np.linspace(0, np.linalg.norm(G_2), dimensions[1])
X, Y = np.meshgrid(x, y)
fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, data2d)
plt.xlabel('G_3')
plt.ylabel('G_2')
fig2 = plt.figure(1)
ay = fig2.add_subplot(111)
ay.set_xlim([0, dimensions[1]/2])
ay.plot(data2d[:, 0])
#plt.show()

#Create 2D array of indices
#indices = list((row_ind,col_ind,K_z) for row_ind in xrange(dimensions[1]) for col_ind in xrange(dimensions[0]))
indices = list((float(K_x)/(dimensions[0]-1)*np.array(G_1) + float(row_ind)/(dimensions[1]-1)*np.array(G_2)+float(col_ind)*np.array(G_1)) for row_ind in xrange(dimensions[1]) for col_ind in xrange(dimensions[0]))
data_diag = np.diag(data2d)
fig3 = plt.figure(2)
az = fig3.add_subplot(111)
az.set_xlim([0, dimensions[1]/2])
az.plot(data_diag)
#plt.show()

#Get corresponding index value for a point in the data grid
#for index,value in ndenumerate(data2d):
#    print(28*index[0]+index[1])

#Write list of data to file for Mathematica
#f = open('data.dat', 'w')
#a = ndenumerate(data2)
#for (i, j, k), v in a:
    #print>>f, i, j, k, v
#f.close

f = open('data_' + filename + '.dat', 'w')
kpointdata = np.zeros(4)
kpoint = np.zeros(3)
finalvec = np.zeros(3)
smallcount = 0
maxZ = 0
minZ = 0

#for i in range(-1*dimensions[0] + 1, dimensions[0] -1, 1):
#    for j in range(-1*dimensions[1] + 1, dimensions[1] - 1, 1):
#        for k in range(-1*dimensions[2] + 1, dimensions[2] - 1, 1):
for i, j, k in itertools.product(xrange(dimensions[0]), xrange(dimensions[1]), xrange(dimensions[2])):
    #print i, j, k
    kpoint = [0, 0, 0]
    kpoint = i*np.array(G_1)/(dimensions[0] - 1) + j*np.array(G_2)/(dimensions[1] - 1) + k*np.array(G_3)/(dimensions[2] - 1)
    normk = np.linalg.norm(kpoint)
    finalvec = [0, 0, 0]
    pluszerominus = [-1, 0, 1]
    for a in pluszerominus:
        for b in pluszerominus:
            for c in pluszerominus:
                testvec = a*np.array(G_1) + b*np.array(G_2) + c*np.array(G_3)
                newvec = kpoint-testvec
                #print(normk, np.linalg.norm(newvec, 1))
                if np.linalg.norm(newvec) < normk:
                    smallcount = smallcount + 1
                    normk = np.linalg.norm(newvec)
                    finalvec = testvec.copy()
                    #print(kpoint, finalvec)
    kpoint = kpoint - finalvec
    kpointdata[0:3] = kpoint
    kpointdata[-1] = data2[i][j][k]
#find the max and min k_z values for plotting multiple BZs
    for test in [-1, 1]:
        if test*kpointdata[2] > maxZ:
            maxZ = test*kpointdata[2]
        if test*kpointdata[2] < minZ:
            minZ = test*kpointdata[2]
    #outstr1 = "{}\t{}\t{}\t{}\n".format(*kpointdata)
    outstr1 = "{}\t{}\t{}\t{}\n".format('%.6f' % (kpointdata[0] + 0.0), '%.6f' % (kpointdata[1] + 0.0), '%.6f' % (kpointdata[2] + 0.0), '%.6f' % (kpointdata[-1] + 0.0))
    '''outstr2 = "{}\t{}\t{}\t{}\n".format(kpointdata[0], kpointdata[1], -1*kpointdata[2] + 0, kpointdata[-1])
    outstr3 = "{}\t{}\t{}\t{}\n".format(kpointdata[0], -1*kpointdata[1] + 0, kpointdata[2] + 0, kpointdata[-1])
    outstr4 = "{}\t{}\t{}\t{}\n".format(kpointdata[0], -1*kpointdata[1] + 0, -1*kpointdata[2] + 0, kpointdata[-1])
    outstr5 = "{}\t{}\t{}\t{}\n".format(-1*kpointdata[0] + 0, -1*kpointdata[1] + 0, kpointdata[2] + 0, kpointdata[-1])
    outstr6 = "{}\t{}\t{}\t{}\n".format(-1*kpointdata[0] + 0, -1*kpointdata[1] + 0, -1*kpointdata[2] + 0, kpointdata[-1])
    outstr7 = "{}\t{}\t{}\t{}\n".format(-1*kpointdata[0] + 0, kpointdata[1], kpointdata[2] + 0, kpointdata[-1])
    outstr8 = "{}\t{}\t{}\t{}\n".format(-1*kpointdata[0] + 0, kpointdata[1], -1*kpointdata[2] + 0, kpointdata[-1])'''
    #f.write(outstr1 + outstr2 + outstr3 + outstr4 + outstr5 + outstr6 + outstr7 + outstr8)
    f.write(outstr1)
f.close()
#print(maxZ, minZ)
