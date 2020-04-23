# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import pyfea.fea as fea
from idealab_tools.data_exchange import dat
import numpy
import os


def u_d(x):
    mm = x.shape[0]
    M = numpy.zeros((3*mm,3))
    W = numpy.zeros((3*mm,1))
    M[0::3,0] = 1
    M[1::3,1] = 1
    M[2::3,2] = 1
    aa = x[:,1]<-50
    bb = aa.nonzero()[0]
    W[3*bb+2] = .1
    return W,M

material = fea.Material(100000,.3)
factor=100


if __name__=='__main__':
    directory = './'
else:
    import sys
    m = sys.modules[__name__]
    directory = os.path.split(m.__file__)[0]


coordinates = dat.read(os.path.join(directory,'coordinates.dat'),float)
elements = dat.read(os.path.join(directory,'elements.dat'),int) - 1
dirichlet = dat.read(os.path.join(directory,'dirichlet.dat'),int) - 1
neumann = dat.read(os.path.join(directory,'neumann.dat'),int) - 1
tris = numpy.r_[dirichlet,neumann]

dirichlet_nodes = numpy.unique(dirichlet)
neumann_nodes = numpy.unique(neumann)

fea.plot_triangles(coordinates,tris)

x,u = fea.compute(material,coordinates,elements,[],neumann,dirichlet_nodes,fea.volume_force_empty,fea.surface_force_empty,u_d)
ax = fea.show(elements,[],tris,coordinates,u,material,factor=factor) 
fea.plot_nodes(coordinates,dirichlet_nodes,u,ax,factor)

output= {}
output['x']=x

import pyfea.error_check as error_check
import idealab_tools.data_exchange.dat

#filename = os.path.join(directory,'bracket.yaml')
#error_check.error_check(output,filename,generate_new = False)        

x2 = idealab_tools.data_exchange.dat.read('results/x.dat',float)

for key,value in output.items():
    dat_filename = 'results/'+key+'.dat'
    a = error_check.compare_matrices(value,dat_filename,tol=1e-6)
    if a>0:
        raise(Exception('too many errors'))

#import idealab_tools.plot_tris as pt
#pt.plot_tris(xyz,triangles,verts_colors = colors, draw_edges = True, edge_color=(0,0,0,1))


