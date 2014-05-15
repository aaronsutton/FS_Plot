#import re
import numpy as np
import matplotlib.pyplot as plt
import itertools
#import matplotlib.figure as fig
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
# from tvtk.api import tvtk, write_data
# from scipy.interpolate import griddata
from pyhull.voronoi import VoronoiTess
from return_data import return_data
#import gts
#from pyhull.delaunay import DelaunayTri
#import mayavi
#import pyvtk

assert plt
assert Axes3D
assert itertools
# assert tvtk
#assert pyvtk

filename = 'YbSO25072013.bxsf.band-39'
#Try maya plotting stuff
#f_plot = open('data_' + filename + '.dat', 'r')
#plotData = np.loadtxt(f_plot)
#plotX = plotData[:, 0]
#plotY = plotData[:, 1]
#plotZ = plotData[:, 2]
#Energy = plotData[:, 3]
#R = plotData[:, 0:3]
#pts = 30j
#X, Y, Z = np.mgrid[min(plotX):max(plotX):pts, min(plotY):max(plotY):pts, min(plotZ):max(plotZ):pts]
#F = griddata(R, Energy, (X, Y, Z), method='linear')
#f_plot.close()
X, Y, Z, F = return_data(filename)

NN_plot = open('NNgen.dat', 'r')
points = np.loadtxt(NN_plot)
#NN_plot.close()
pointsarray = np.array(points)
#points3d(pointsarray[:, 0], pointsarray[:, 1], pointsarray[:, 2], color=(1, 0, 0))
#points3d(0, 0, 0, color=(1, 0, 0))
d = VoronoiTess(points)
dpoints = d.points
ridges = d.ridges.items()
vertices = d.vertices
vertarray = np.array(vertices)
maxvert = max(vertarray[:, 2])
bzfinalx = list()
bzfinaly = list()
bzfinalz = list()
bzlines = open('bzlines.dat', 'w')
#mayawindow = mayavi.mayavi()i

mlab.figure(1, bgcolor=(1,1,1), size=(325, 325))

for (p1, p2), vertind in ridges:
    if p1 == 0:
        bzplotx = list()
        bzploty = list()
        bzplotz = list()
        for aind in vertind:
            bzplotx.append(vertices[aind][0])
            bzploty.append(vertices[aind][1])
            bzplotz.append(vertices[aind][2])
        #This was used to create closed faces for the outline of the BZ, not needed if creating a Delaunay 3D surface
        if bzplotx[0] != bzplotx[-1] or bzploty[0] != bzploty[-1] or bzplotz[0] != bzplotz[-1]:
            bzplotx.append(bzplotx[0])
            bzploty.append(bzploty[0])
            bzplotz.append(bzplotz[0])
        bzfinalx = bzfinalx + bzplotx
        bzfinaly = bzfinaly + bzploty
        bzfinalz = bzfinalz + bzplotz
        bzplotxa = np.array(bzplotx) * 1
        bzplotya = np.array(bzploty) * 1
        bzplotza = np.array(bzplotz) * 1
        '''for i in range(len(bzplotxa)):
            outline = "{}\t{}\t{}\n".format(bzplotxa[i], bzplotya[i], bzplotza[i])
            bzlines.write(outline)'''
        #Plot the outline
        mlab.plot3d(bzplotxa, bzplotya, bzplotza, color=(0, 0, 0), line_width=0.001, tube_radius=0.001)
bzlines.close()
scalefactor = max(bzfinalz)/Z.max()
Xnew = X.copy()
Ynew = Y.copy()
Znew = Z.copy()
#print Znew.max()
#contour3d(X, Y, Z + 2*max(plotZ), F)

bzfile = open('bzpoints.dat', 'w')
for i in range(len(bzfinalx)):
    outstring = "{}\t{}\t{}\n".format(bzfinalx[i], bzfinaly[i], bzfinalz[i])
    bzfile.write(outstring)
bzfile.close()
FS = mlab.contour3d(Xnew, Ynew, Znew, F, contours=[0.765467])
bzpoints = mlab.points3d(bzfinalx, bzfinaly, bzfinalz, color=(0, 0, 0), scale_factor=0.001)

# Plot a second band?
z_shift = 1
if z_shift:
    z_shift_val = Z.max() - Z.min()
    Z_z = [x + z_shift_val for x in Z]
    FS_Z = mlab.contour3d(X, Y, Z_z, F, contours=[0.765467])
    xz_shift_val = 2*X.max()
    z_shift_val = Z.max()
    X_xz = [x + xz_shift_val - 0.02 for x in X]
    Z_xz = [x + z_shift_val for x in Z]
    FS_XZ = mlab.contour3d(X_xz, Y, Z_xz, F, contours=[0.765467])

plot_band2 = 0
filename2 = 'YbSO25072013.bxsf.band-40'

if plot_band2:
    X, Y, Z, F = return_data(filename2)
    FS2 = mlab.contour3d(X, Y, Z, F, contours=[0.765467])

# Save plots - various options
#mlab.savefig(filename='test.wrl')
#savefig('FS' + filename + '.png')
#delmesh = mlab.pipeline.delaunay3d(bzpoints)
#bzpoints.remove
#delsurface = mlab.pipeline.surface(delmesh, opacity=0.5)




# This was a test of using the gts package. It sort of works, but isn't necessary
'''isotest = gts.isosurface(F,0.765467)
x,y,z,t = gts.get_coords_and_face_indices(isotest,True)
mlab.triangular_mesh(x,y,z,t)
v = list()
for [a,b,c] in vertices:
    v.append(gts.Vertex(a,b,c))
outerlist = list()
for (p1, p2), vertind in ridges:
    if p1 == 0:
        innerlist = list()
        for countvert in range(len(vertind)):
            if countvert < len(vertind) - 1:
                innerlist.append(gts.Edge(v[vertind[countvert]],v[vertind[countvert + 1]]))
            else:
                innerlist.append(gts.Edge(v[vertind[countvert]],v[vertind[0]]))
        outerlist.append(innerlist)
faces = list()
for o in outerlist:
    facesinner = list()
    for p in o:
        facesinner.append(p)
    print facesinner'''




# clip = mlab.pipeline.data_set_clipper(FS)
# transform = tvtk.Transform()
# transform.scale(0.9, 1, 1.2)
# surf = mlab.pipeline.iso_surface(clip)
# transform.matrix.__setstate__({'elements':
#         [ 0.9,  0.1,  0.1, -1.3,
#           0.0,  0.9, -0.4,  2.9,
#          -0.1,  0.3,  0.9, -0.6,
#           0.0,  0.0,  0.0,  1.0]
#      })
# clip.widget.widget.set_transform(transform)
# clip.widget.update_implicit_function()
