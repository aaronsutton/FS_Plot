import numpy as np
#import matplotlib.pyplot as plt
#import itertools

with open('/home/asutton/Dropbox/WIEN2K/CASES/YbSO25072013/YbSO25072013.bxsf.band-39', 'r') as f_id:
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
    #print 'Number of bands is {}'.format(no_bands)
    line = f_id.readline()
    dimensions = int(line.split()[0]), int(line.split()[1]), int(line.split()[2])
    line = f_id.readline()
    line = f_id.readline()
    G_1 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
    line = f_id.readline()
    G_2 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
    line = f_id.readline()
    G_3 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])


f = open('NNgen.dat', 'w')
#a = 1.0/0.13196400
#c = 1.0/0.05377300
#a = 4.01000
#c = 9.84100
#G_1 = [1.0/a, 0, 1.0/c]
#G_2 = [0, 1.0/a, 1.0/c]
#G_3 = [1.0/a, 1.0/a, 0]
#print(G_1)
i, j, k = 0, 0, 0
plusminusone = [0, -1, 1]
for i in plusminusone:
            for j in plusminusone:
                for k in plusminusone:
                	NN = i*np.array(G_1) + j*np.array(G_2) + k*np.array(G_3)
                	outstr = "{}\t{}\t{}\n".format(NN[0], NN[1], NN[2])
                	f.write(outstr)
f.close()
