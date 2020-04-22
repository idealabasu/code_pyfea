# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import numpy
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import idealab_tools.fea_tetra.fea as fea

x = numpy.r_[0:1:10j]
y = numpy.r_[0:.1:3j]
z = numpy.r_[0:.1:3j]
xyz = numpy.array(numpy.meshgrid(x,y,z)).T
xyz = xyz.reshape((-1,3))

d = Delaunay(xyz)
xyz = d.points
quad = d.simplices
quad = fea.fix_tet_order(xyz,quad)
quad = fea.remove_zero_volume(xyz,quad)

tri_indeces = numpy.c_[[0,1,2],[1,2,3],[2,3,0],[3,0,1]]
tris = quad[:,tri_indeces]
tris = tris.transpose((0,2,1)).reshape((-1,3))
tris = numpy.array(list(set([tuple(sorted(item)) for item in tris])))

#fea.plot_triangles(xyz,tris)

xx = xyz[:,0]
yy = xyz[:,1]
zz = xyz[:,2]


material = fea.Material(100000,.3)
factor = 100
coordinates = xyz
elements = quad

z_max = coordinates.max(0)[2]
z_min = coordinates.min(0)[2]
x_max = coordinates.max(0)[0]
x_min = coordinates.min(0)[0]


ii_bottom = ((coordinates[tris,0]==x_min).sum(1)==3)
ii_top = ((coordinates[tris,0]==x_max).sum(1)==3)
dirichlet_bottom = tris[ii_bottom]
dirichlet_top = tris[ii_top]
dirichlet = numpy.r_[dirichlet_bottom,dirichlet_top]

neumann = numpy.zeros((0,3),dtype = int)
dirichlet_nodes = numpy.unique(dirichlet)
neumann_nodes = numpy.unique(neumann)

#def u_d(x):
#    mm = x.shape[0]
#    M = numpy.zeros((3*mm,3))
#    W = numpy.zeros((3*mm,1))
#    M[0::3,0] = 1
#    M[1::3,1] = 1
#    M[2::3,2] = 1
#    aa = x[:,0]==1
#    bb = aa.nonzero()[0]
#    W[3*bb+2] = 1e-3
#    return W,M

def u_d(x):
    mm = x.shape[0]
    M = numpy.zeros((3*mm,3))
    W = numpy.zeros((3*mm,1))
    
    aa = (x[:,0]==1).nonzero()[0]
    bb = (x[:,0]!=1).nonzero()[0]

    M[3*bb,0] = 1
    M[3*bb+1,1] = 1
    M[3*bb+2,2] = 1

    M[3*aa+2,2] = 1
    W[3*aa+2] = 1e-3
    return W,M

#
x,u = fea.compute(material,coordinates,elements,[],neumann,dirichlet_nodes,fea.volume_force_empty,fea.surface_force_empty,u_d)
ax = fea.show(elements,[],tris,coordinates,u,material,factor=factor) 
fea.plot_nodes(coordinates,dirichlet_nodes,u,ax,factor)

filename = 'cantilever_delaunay.yaml'

output = {}
output['x']=x
output['u']=u

import idealab_tools.fea_tetra.error_check as error_check
error_check.error_check(output,filename,generate_new = False)