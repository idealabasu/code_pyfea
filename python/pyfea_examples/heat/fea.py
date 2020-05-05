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

def VolumeForce(x,t):
    VolumeForce = numpy.zeros((1,1))
    return VolumeForce

def Stress(x,t):
    Stress = numpy.zeros((1,1))
    return Stress

def stima3(vertices):
    d = vertices.shape[1]
    
    A = numpy.r_[numpy.ones((1,d+1)),vertices.T]
    b = numpy.r_[numpy.zeros((1,d)),numpy.eye(d)]
    G = scipy.linalg.solve(A,b)
    C = scipy.linalg.det(A)
    M = C * G.dot(G.T) / (numpy.r_[1:d+1].prod())

    return M

def DirichletBoundaryValue(x,t):
    f1 = x[:,0]==0
    f2 = x[:,1]==3
    DirichletBoundaryValue =  numpy.zeros((x.shape[0],1))
    DirichletBoundaryValue[f1] = 1
    DirichletBoundaryValue[f2] = 1
    return DirichletBoundaryValue

coordinates = dat.read(os.path.join(directory,'coordinates.dat'),float)
coordinates = coordinates[:,1:]
elements3 = dat.read(os.path.join(directory,'elements3.dat'),float)
elements3 = numpy.array(elements3,dtype=int)
elements3 = elements3[:,1:]
elements3 -=1
dirichlet = dat.read(os.path.join(directory,'dirichlet.dat'),float)
dirichlet = numpy.array(dirichlet,dtype=int)
dirichlet = dirichlet[:,1:]
dirichlet -=1
neumann = dat.read(os.path.join(directory,'neumann.dat'),float) - 1
neumann = numpy.array(neumann ,dtype=int)
neumann = neumann[:,1:]
neumann -=1

m = coordinates.shape[0]

free_nodes=sorted(set(range(m))-set(dirichlet.flatten().tolist()))

A = numpy.zeros((m,m))
B = numpy.zeros((m,m))

T = 1
dt = 0.01
N = int(T/dt)

U = numpy.zeros((m,N+1))
U0 = 0

for item in elements3:
    kk,ll = numpy.c_[numpy.meshgrid(item,item)].transpose(0,2,1).reshape(2,-1)
    A[kk,ll] += stima3(coordinates[item,:]).flatten()
    B[kk,ll] += (scipy.linalg.det(numpy.r_[numpy.array([[1,1,1]]),coordinates[item,:].T])*(numpy.array([[2,1,1],[1,2,1],[1,1,2]])/24)).flatten()
    
U[:,0] = U0

dirichlet_nodes = numpy.unique(dirichlet)

for nn in numpy.r_[1:N+1]:
    b = numpy.zeros((m,1))
    for item in elements3:
        bb = scipy.linalg.det(numpy.r_[numpy.array([[1,1,1]]),coordinates[item,:].T]) * dt*VolumeForce(coordinates[item,:].sum(0)/3,nn*dt)/6
        b[item] += bb

    for item in neumann:
        bb = scipy.linalg.norm(coordinates[item[0],:]-coordinates[item[1],:]) * dt*Stress((coordinates[item,:]).sum(0)/2,nn*dt)/2
        b[item] += bb

    b += (B.dot(U[:,nn-1]))[:,None]

    u = numpy.zeros((m,1))
    u[dirichlet_nodes] = DirichletBoundaryValue(coordinates[dirichlet_nodes,:],nn*dt)
    
    b -= (dt*A + B).dot(u)

    oo = len(free_nodes)
    AA = numpy.zeros((oo,oo))
    kk,ll = numpy.c_[numpy.meshgrid(free_nodes,free_nodes)].transpose(0,2,1).reshape(2,-1)
    AAA = (dt*A[kk,ll] + B[kk,ll])
    kk,ll = numpy.c_[numpy.meshgrid(range(oo),range(oo))].transpose(0,2,1).reshape(2,-1)
    AA[kk,ll] = AAA
    bb = b[free_nodes]
    u[free_nodes] = scipy.linalg.solve(AA,bb)
    U[:,nn] = u.flatten()

output= {}
output['U']=U

import pyfea.error_check as error_check
import idealab_tools.data_exchange.dat

for key,value in output.items():
    dat_filename = 'results/'+key+'.dat'
    a = error_check.compare_matrices(value,dat_filename)
    if a>0:
        raise(Exception('too many errors'))

uu = U[:,-1]
coordinates3 = numpy.c_[coordinates,uu]
vc = numpy.array([cm.jet(item) for item in uu])
#idealab_tools.plot_tris.plot_tris(coordinates3,elements3,face_colors = (1,0,0,1), draw_edges = True,edge_color=(0,0,0,1))
idealab_tools.plot_tris.plot_tris(coordinates3,elements3,verts_colors = vc, draw_edges = True,edge_color=(0,0,0,1))

