# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 13:38:50 2018

@author: danaukes
"""
import numpy
import scipy.linalg
import os
from idealab_tools.data_exchange import dat
import idealab_tools.plot_tris
import matplotlib.cm as cm
directory = ''
import pyfea.fea as fea
fea.NDIM=2
fea.NBOUNDARYNODES=2

#def VolumeForce(x,t):
#    VolumeForce = numpy.zeros((1,1))
#    return VolumeForce
#
#def Stress(x,t):
#    Stress = numpy.zeros((1,1))
#    return Stress
#
#def stima3(vertices):
#    d = vertices.shape[1]
#    
#    A = numpy.r_[numpy.ones((1,d+1)),vertices.T]
#    b = numpy.r_[numpy.zeros((1,d)),numpy.eye(d)]
#    G = scipy.linalg.solve(A,b)
#    C = scipy.linalg.det(A)
#    M = C * G.dot(G.T) / (numpy.r_[1:d+1].prod())
#
#    return M
#
#def DirichletBoundaryValue(x,t):
#    f1 = x[:,0]==0
#    f2 = x[:,1]==3
#    DirichletBoundaryValue =  numpy.zeros((x.shape[0],1))
#    DirichletBoundaryValue[f1] = 1
#    DirichletBoundaryValue[f2] = 1
#    return DirichletBoundaryValue


def localdj(vertices,U):
    U = U.flatten()
    Eps = 1/100
    gdem = numpy.concatenate([numpy.ones((1,3)),vertices.T],0)
    gnum = numpy.concatenate([numpy.zeros((1,2)),numpy.eye(2)],0)
    G = numpy.linalg.solve(gdem,gnum)
    Area = numpy.linalg.det(gdem)/ 2;
    fac1 = numpy.array([[2,1,1],[1,2,1],[1,1,2]])/12 
    fac2 = numpy.array([[12*U[1-1]**2+2*(U[2-1]**2+U[3-1]**2+U[2-1]*U[3-1])+6*U[1-1]*(U[2-1]+U[3-1]),
      3*(U[1-1]**2+U[2-1]**2)+U[3-1]**2+4*U[1-1]*U[2-1]+2*U[3-1]*(U[1-1]+U[2-1]),
      3*(U[1-1]**2+U[3-1]**2)+U[2-1]**2+4*U[1-1]*U[3-1]+2*U[2-1]*(U[1-1]+U[3-1])],
  [3*(U[1-1]**2+U[2-1]**2)+U[3-1]**2+4*U[1-1]*U[2-1]+2*U[3-1]*(U[1-1]+U[2-1]),
      12*U[2-1]**2+2*(U[1-1]**2+U[3-1]**2+U[1-1]*U[3-1])+6*U[2-1]*(U[1-1]+U[3-1]),
      3*(U[2-1]**2+U[3-1]**2)+U[1-1]**2+4*U[2-1]*U[3-1]+2*U[1-1]*(U[2-1]+U[3-1])],
  [3*(U[1-1]**2+U[3-1]**2)+U[2-1]**2+4*U[1-1]*U[3-1]+2*U[2-1]*(U[1-1]+U[3-1]),
      3*(U[2-1]**2+U[3-1]**2)+U[1-1]**2+4*U[2-1]*U[3-1]+2*U[1-1]*(U[2-1]+U[3-1]),
      12*U[3-1]**2+2*(U[1-1]**2+U[2-1]**2+U[1-1]*U[2-1])+6*U[3-1]*(U[1-1]+U[2-1])]])/60
    M = Area*(Eps*G.dot(G.T)-fac1+fac2)
    return M
    
def localj(vertices,U):
    U = U.flatten()
    Eps = 1/100;
    gdem = numpy.concatenate([numpy.ones((1,3)),vertices.T],0)
    gnum = numpy.concatenate([numpy.zeros((1,2)),numpy.eye(2)],0)
    G = numpy.linalg.solve(gdem,gnum)
    Area = numpy.linalg.det(gdem)/ 2;
    fac1 = numpy.array([[2,1,1],[1,2,1],[1,1,2]])/12 
    fac2 = numpy.array([[4*U[1-1]**3+ U[2-1]**3+U[3-1]**3+3*U[1-1]**2*(U[2-1]+U[3-1])+2*U[1-1]*(U[2-1]**2+U[3-1]**2)+U[2-1]*U[3-1]*(U[2-1]+U[3-1])+2*U[1-1]*U[2-1]*U[3-1]],
      [4*U[2-1]**3+ U[1-1]**3+U[3-1]**3+3*U[2-1]**2*(U[1-1]+U[3-1])+2*U[2-1]*(U[1-1]**2+U[3-1]**2)+U[1-1]*U[3-1]*(U[1-1]+U[3-1])+2*U[1-1]*U[2-1]*U[3-1]],
      [4*U[3-1]**3+ U[2-1]**3+U[1-1]**3+3*U[3-1]**2*(U[2-1]+U[1-1])+2*U[3-1]*(U[2-1]**2+U[1-1]**2)+U[2-1]*U[1-1]*(U[2-1]+U[1-1])+2*U[1-1]*U[2-1]*U[3-1]]])/60
    fac2 = fac2.flatten()
    b=Area*((Eps*G.dot(G.T)- fac1).dot(U)+ fac2)
    return b 

def f(x):
    volume_force = numpy.zeros((x.shape[0],1))
    return volume_force

def g(x):
    stress = numpy.zeros((x.shape[0],1))
    return stress

coordinates = dat.read(os.path.join(directory,'coordinates.dat'),float)
coordinates = coordinates[:,1:]
elements3 = numpy.array(dat.read(os.path.join(directory,'elements3.dat'),float),dtype=numpy.int) - 1
elements3 = elements3[:,1:]
dirichlet = numpy.array(dat.read(os.path.join(directory,'dirichlet.dat'),float),dtype = numpy.int) - 1
dirichlet = dirichlet[:,1:]
#neumann = numpy.array(dat.read(os.path.join(directory,'neumann.dat'),float),dtype = numpy.int) - 1
#neumann = neumann[:,1:]
neumann= numpy.zeros((0,2))

m = coordinates.shape[0]

free_nodes=sorted(set(range(m))-set(dirichlet.flatten().tolist()))

U = -numpy.ones((m,1))

for ii in range(50):
    A = numpy.zeros((m,m))
    b = numpy.zeros((m,1))
    
    for item in elements3:
        coords = coordinates[item,:]
        
        kk,ll = numpy.c_[numpy.meshgrid(item,item)].transpose(0,2,1).reshape(2,-1)
        ldj = localdj(coords,U[item])
        A[kk,ll] += ldj.flatten()
        
        lj = numpy.array([localj(coords,U[item])]).T
        b[item] += lj
        
        aa =  numpy.concatenate([[[1,1,1]],coords.T],0).T
        bb = numpy.linalg.det(aa)
        cc = bb*f(numpy.array([coords.sum(0)])/3)/6
        b[item] += cc
    
    for item in neumann:
        coords = coordinates[item,:]
        dd = (numpy.linalg.norm(coords) + coords.dot(g(numpy.array([coords.sum(0)])/2)/2))
        b[item] -= dd


    #Dirichlet conditions
    W = numpy.zeros((coordinates.shape[0],1))
    W[numpy.unique(dirichlet)] = 0
    
    #Solving one Newton step
    n = len(free_nodes)
    kk,ll = numpy.c_[numpy.meshgrid(free_nodes,free_nodes)].transpose(0,2,1).reshape(2,-1)
    A_free = A[kk,ll]
    A_free = A_free.reshape(n,n)
    b_free = b[free_nodes]
    W[free_nodes] = numpy.linalg.solve(A_free,b_free)
    U -= W
    if numpy.linalg.norm(W) < 1e-10:
        break

#material = fea.Material(1e7,.3)

import idealab_tools.plot_tris as pt
pt.plot_tris(coordinates,elements3,face_colors = (1,0,0,1), draw_edges = True, edge_color=(0,0,0,1))

#

#    T = 1
#    dt = 0.01
#    N = int(T/dt)
#    
#    U0 = 0
#    
#    for item in elements3:
#        kk,ll = numpy.c_[numpy.meshgrid(item,item)].transpose(0,2,1).reshape(2,-1)
#        A[kk,ll] += stima3(coordinates[item,:]).flatten()
#        B[kk,ll] += (scipy.linalg.det(numpy.r_[numpy.array([[1,1,1]]),coordinates[item,:].T])*(numpy.array([[2,1,1],[1,2,1],[1,1,2]])/24)).flatten()
#        
#    U[:,0] = U0
#    
#    dirichlet_nodes = numpy.unique(dirichlet)
#    
#    for nn in numpy.r_[1:N+1]:
#        b = numpy.zeros((m,1))
#        for item in elements3:
#            bb = scipy.linalg.det(numpy.r_[numpy.array([[1,1,1]]),coordinates[item,:].T]) * dt*VolumeForce(coordinates[item,:].sum(0)/3,nn*dt)/6
#            b[item] += bb
#    
#        for item in neumann:
#            bb = scipy.linalg.norm(coordinates[item[0],:]-coordinates[item[1],:]) * dt*Stress((coordinates[item,:]).sum(0)/2,nn*dt)/2
#            b[item] += bb
#    
#        b += (B.dot(U[:,nn-1]))[:,None]
#    
#        u = numpy.zeros((m,1))
#        u[dirichlet_nodes] = DirichletBoundaryValue(coordinates[dirichlet_nodes,:],nn*dt)
#        
#        b -= (dt*A + B).dot(u)
#    
#        oo = len(free_nodes)
#        AA = numpy.zeros((oo,oo))
#        kk,ll = numpy.c_[numpy.meshgrid(free_nodes,free_nodes)].transpose(0,2,1).reshape(2,-1)
#        AAA = (dt*A[kk,ll] + B[kk,ll])
#        kk,ll = numpy.c_[numpy.meshgrid(range(oo),range(oo))].transpose(0,2,1).reshape(2,-1)
#        AA[kk,ll] = AAA
#        bb = b[free_nodes]
#        u[free_nodes] = scipy.linalg.solve(AA,bb)
#        U[:,nn] = u.flatten()
#    
#    output= {}
#    output['U']=U
#    
#    import pyfea.error_check as error_check
#    import idealab_tools.data_exchange.dat
#    
#    #for key,value in output.items():
#    #    dat_filename = 'results/'+key+'.dat'
#    #    a = error_check.compare_matrices(value,dat_filename)
#    #    if a>0:
#    #        raise(Exception('too many errors'))
#    error_check.compare_output(output,'results')
#    
#    uu = U[:,-1]
#    coordinates3 = numpy.c_[coordinates,uu]
#    vc = numpy.array([cm.jet(item) for item in uu])
#    #idealab_tools.plot_tris.plot_tris(coordinates3,elements3,face_colors = (1,0,0,1), draw_edges = True,edge_color=(0,0,0,1))
#    idealab_tools.plot_tris.plot_tris(coordinates3,elements3,verts_colors = vc, draw_edges = True,edge_color=(0,0,0,1))
#    
