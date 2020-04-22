# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import idealab_tools.fea_tetra.fea as fea
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
    directory = './bracket'
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

x,u = fea.compute(material,coordinates,elements,neumann,dirichlet_nodes,fea.volume_force_empty,fea.surface_force_empty,u_d)
ax = fea.show(elements,tris,coordinates,u,material,factor=factor) 
fea.plot_nodes(coordinates,dirichlet_nodes,u,ax,factor)


filename = os.path.join(directory,'bracket.yaml')

output = {}
output['x']=x
output['u']=u

import idealab_tools.fea_tetra.error_check as error_check
error_check.error_check(output,filename,generate_new = False)
    




#    errors = 0 
#    import idealab_tools.fea_tetra.error_check as error_check
#    errors += error_check.compare_matrices(b,'b.dat','results')
#    errors += error_check.compare_matrices(B,'B_big.dat','results')
#    errors += error_check.compare_matrices(W,'W.dat','results')
#    errors += error_check.compare_matrices(M,'M.dat','results')
#    errors += error_check.compare_matrices(u,'u.dat','results')
#    errors += error_check.compare_matrices(x,'x.dat','results')
#    errors += error_check.compare_matrices(A,'A.dat','results')
#    print('num errors: '+str(errors))

