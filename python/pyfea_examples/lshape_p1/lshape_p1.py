# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
-"""

import pyfea.fea as fea
fea.NDIM=2
fea.NBOUNDARYNODES=2
from idealab_tools.data_exchange import dat
import numpy
import os
from math import pi
from numpy import sin,cos,arctan2

material = fea.Material(100000,.3)
factor=100

def cart2pol_2d(x,y):
    r = (x**2+y**2)**.5
    theta = arctan2(y,x)
    return r,theta

def u_value(x,Lambda,mu):
    r,phi = cart2pol_2d(x[:,0],x[:,1])
    alph = .544483737
    omega = 3*pi/4
    C_1 = -cos((alph+1)*omega)/cos((alph-1)*omega)
    C_2 = 2*(Lambda+2*mu)/(Lambda+mu)
    ut = (1/(2*mu)) * r**alph *((alph+1)*sin((alph+1)*phi)+ (C_2+alph-1)*C_1*sin((alph-1)*phi))
    ur = (1/(2*mu))*r**alph * (-(alph+1)*cos((alph+1)*phi)+ (C_2-(alph+1))*C_1*cos((alph-1)*phi))
    value = ur * cos(phi) - ut * sin(phi), ur * sin(phi) + ut * cos(phi)
    return value

def u_d(x):
    mm = x.shape[0]
    M = numpy.zeros((fea.NDIM*mm,fea.NDIM))
    W = numpy.zeros((fea.NDIM*mm,1))
    
    value = u_value(x,material.Lambda,material.mu)

    M[::2,0] = 1
    M[1::2,1] = 1

    W[::2,0] = value[0]
    W[1::2,0] = value[1]
    
    return W,M

directory = './'

def surface_force(x,n):
    sforce = numpy.zeros((x.shape[0],fea.NDIM))
    ii = (n[:,1]==1).nonzero()[0]
    sforce[ii,1]=1
    return sforce

coordinates = dat.read(os.path.join(directory,'coordinates.dat'),float)
elements3 = dat.read(os.path.join(directory,'elements3.dat'),int) - 1
elements4 = numpy.zeros((0,4))
dirichlet = dat.read(os.path.join(directory,'dirichlet.dat'),int) - 1
neumann = numpy.zeros((0,2))
tris = numpy.r_[dirichlet,neumann]

dirichlet_nodes = numpy.unique(dirichlet)
neumann_nodes = numpy.unique(neumann)

#fea.plot_triangles(coordinates,tris)

x,u = fea.compute(material,coordinates,elements3,elements4,neumann,dirichlet_nodes,fea.volume_force_empty,surface_force,u_d)

Eps3 = numpy.zeros((elements3.shape[0],4))
Sigma3 = numpy.zeros((elements3.shape[0],4))
Eps4 = numpy.zeros((elements4.shape[0],4))
Sigma4 = numpy.zeros((elements4.shape[0],4))
AreaOmega = numpy.zeros((coordinates.shape[0],1))
AvS = numpy.zeros((coordinates.shape[0],4))
AvE = numpy.zeros((coordinates.shape[0],4))






for jj,row in enumerate(elements3):
    augmented = fea.augment_2d(coordinates[row])
    area3 = numpy.linalg.det(augmented)/2
    AreaOmega[row] += area3

    PhiGrad = fea.phigrad(augmented)
    ii = 2*numpy.array([[1,1]]).T.dot(numpy.array([row])) + numpy.array([[0,1]]).T.dot(numpy.array([[1,1,1]]))
    ui = u[ii].squeeze()
    U_Grad = ui.dot(PhiGrad )

    Eps3[jj] = ((U_Grad+U_Grad.T)/2).reshape(1,4)
    Sigma3[jj] = (material.Lambda*(U_Grad.trace())*numpy.eye(2)+2*material.mu*(U_Grad+U_Grad.T)/2).reshape((1,4))
    AvE[row] += area3*numpy.array([[1,1,1]]).T.dot(Eps3[[jj]])
    AvS[row] += area3*numpy.array([[1,1,1]]).T.dot(Sigma3[[jj]])

for jj,row in enumerate(elements4):
    augmented = fea.augment_2d(coordinates[row[:3]])
    area4 = numpy.linalg.det(augmented)
    AreaOmega[row] += area4

    Vertices = coordinates[row]
    DD = numpy.array([Vertices[1]-Vertices[0], Vertices[3]-Vertices[0]])
    AA = numpy.linalg.inv(DD)
    CC = numpy.concatenate([u[2*row],u[2*row+1]],1)
    BB = numpy.array([[-1,1,1,-1],[-1,-1,1,1]])
    U_Grad = AA.dot(BB.dot(CC))/2

    Eps4[jj] = ((U_Grad+U_Grad.T)/2).reshape(1,4)
    Sigma4[jj] = (material.Lambda*(U_Grad.trace())*numpy.eye(2)+2*material.mu*(U_Grad+U_Grad.T)/2).reshape((1,4))
    
    AvE[row] += area4*numpy.array([[1,1,1,1]]).T.dot(Eps4[[jj]])
    AvS[row] += area4*numpy.array([[1,1,1,1]]).T.dot(Sigma4[[jj]])

AvE/=AreaOmega*(numpy.array([1,1,1,1]))
AvS/=AreaOmega*(numpy.array([1,1,1,1]))

s11 = AvS[:,0]
s22 = AvS[:,3]
s12 = AvS[:,1]
mu = material.mu
Lam = material.Lambda
devs2 = (mu**2/(6*(mu+Lam)**2)+1/2)*(s11+s22)**2 + 2*(s12**2-s11*s22)
AvC = numpy.array([devs2/(4*mu)]).T


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.cm

import pyfea.mesh_tools as mt
colors = mt.build_rgba(AvC,matplotlib.cm.gray,invert=True,normalize=True)

triangles = elements3

import pyfea.mesh_tools as mt
xyz = mt.convert_2d_to_3d(coordinates+factor*u.reshape((-1,2)))

output= {}
output['AvS']=AvS
output['AvE']=AvE
#output['u']=u
output['x']=x
output['Sigma3']=Sigma3
#output['Sigma4']=Sigma4
output['Eps3']=Eps3
#output['Eps4']=Eps4
#output['AreaOmega']=AreaOmega
#output['AvC']=AvC

import pyfea.error_check as error_check
import idealab_tools.data_exchange.dat

#filename='lshape_p1.yaml'
#error_check.error_check(output,filename,generate_new = False)        

for key,value in output.items():
    dat_filename = 'results/'+key+'.dat'
    a = error_check.compare_matrices(value,dat_filename)
    if a>0:
        raise(Exception('too many errors'))

#import idealab_tools.plot_tris as pt
#pt.plot_tris(xyz,triangles,verts_colors = colors, draw_edges = True, edge_color=(0,0,0,1))



