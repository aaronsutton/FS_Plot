import numpy as np

def genNN(filename, write_NN):
    with open(filename, 'r') as f_id:
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
        line = f_id.readline()
        dimensions = int(line.split()[0]), int(line.split()[1]), int(line.split()[2])
        line = f_id.readline()
        line = f_id.readline()
        G_1 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
        line = f_id.readline()
        G_2 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
        line = f_id.readline()
        G_3 = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])

    if write_NN:
        f = open('NNgen.dat', 'w')
        i, j, k = 0, 0, 0
        plusminusone = [0, -1, 1]
        for i in plusminusone:
                    for j in plusminusone:
                        for k in plusminusone:
                            NN = i*np.array(G_1) + j*np.array(G_2) + k*np.array(G_3)
                            outstr = "{}\t{}\t{}\n".format(NN[0], NN[1], NN[2])
                            f.write(outstr)
        f.close()
    return E_F, G_1, G_2, G_3
