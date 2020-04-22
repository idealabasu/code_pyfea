# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import pygmsh as pg
import numpy as np
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yaml
from idealab_tools.data_exchange.generic_data import GenericData

import pyfea.fea as fea
import idealab_tools.matplotlib_tools

def select_tri_from_nodes(points,tris,points_filter):
    return points_filter[tris].sum(1)==3

def select_quad_from_nodes(points,quads,points_filter):
    return points_filter[quads].sum(1)==4

def create_t(points,quads):
    quad_points = points[quads,:]
    dqp = quad_points[1:] - quad_points[0]


import pyfea.mesh as pm

meshname = 'plate.yaml'
#mo=load_mesh(meshname)
    

def generate_plate():
    geom = pg.opencascade.Geometry(characteristic_length_max=.002, characteristic_length_min=.001)
    geom.add_box((0,0,0),(.10,.01,.01))
    mo = pg.generate_mesh(geom)
    return mo


mo = pm.get_mesh(meshname,generate_plate,force_regen=False)

triangles_outer = mo.cells['triangle']


material = fea.Material(1e7,.3)
factor = 1

coordinates = mo.points[:]
elements = mo.cells['tetra']

used_elements = fea.find_used_elements(elements,triangles_outer)
coordinates,mapping = fea.coord_reduce(coordinates,used_elements)
triangles_outer = fea.element_reduce(triangles_outer,mapping)
elements= fea.element_reduce(elements,mapping)

a=fea.analyze(coordinates,elements)
print(a)
elements[a] = elements[a][:,(0,2,1,3)]
a=fea.analyze(coordinates,elements)
print(a)

T = coordinates[elements[:,1:]]-coordinates[elements[:,0:1]]
dt = numpy.array([numpy.linalg.det(item) for item in T])
elements = elements[dt!=0]


xx = coordinates[:,0]
yy = coordinates[:,1]
zz = coordinates[:,2]

x_min = coordinates.min(0)[0]
x_max = coordinates.max(0)[0]
y_min = coordinates.min(0)[1]
y_max = coordinates.max(0)[1]
z_min = coordinates.min(0)[2]
z_max = coordinates.max(0)[2]

ii_tri_x_minus = ((coordinates[triangles_outer,0]==x_min).sum(1)==3)
ii_tri_x_plus = ((coordinates[triangles_outer,0]==x_max).sum(1)==3)
ii_tri_y_minus = ((coordinates[triangles_outer,1]==y_min).sum(1)==3)
ii_tri_y_plus = ((coordinates[triangles_outer,1]==y_max).sum(1)==3)
ii_tri_z_minus = ((coordinates[triangles_outer,2]==z_min).sum(1)==3)
ii_tri_z_plus = ((coordinates[triangles_outer,2]==z_max).sum(1)==3)
#ii_neumann = (ii_bottom+ii_top)==0

constrained_tris = triangles_outer[ii_tri_x_minus]
constrained_nodes = numpy.unique(constrained_tris)

surface_forces= numpy.zeros((0,3),dtype = int)
surface_forces = triangles_outer[ii_tri_x_plus]
neumann_nodes = numpy.unique(surface_forces)

#heat_source_nodes = 
#fea.plot_triangles(coordinates,triangles_outer)


def volume_force(x):
    density = 1000
    volforce = numpy.zeros((x.shape[0],3))
#    volforce[:,2] = -9.81*density
#    volforce[:,2] = -9.81
    return volforce    

f = .1
a = .01*.01
p = f/a

def area_force(x,n):
    f = numpy.array([[0,0,-p]])
    return f

point_forces = []
#point_force_tris = triangles_outer[ii_tri_x_plus]
#point_force_nodes = numpy.unique(point_force_tris)
#ll = len(point_force_nodes)
#point_forces = [(ii,[0,0,-.1/ll]) for ii in point_force_nodes]

#def point_force(x):
#    volforce = numpy.zeros((x.shape[0],3))
#    volforce[:,2] = -9.81*density
#    return volforce    

def u_d(x):
    mm = x.shape[0]
    M = numpy.zeros((3*mm,3))
    W = numpy.zeros((3*mm,1))
    
#    aa = (x[:,0]==1).nonzero()[0]
    bb = (x[:,0]==x_min).nonzero()[0]

    M[3*bb,0] = 1
    M[3*bb+1,1] = 1
    M[3*bb+2,2] = 1

#    M[3*aa+2,2] = 1
#    W[3*aa+2] = 1e-1
    return W,M

elements4 = []

x,u = fea.compute(material,coordinates,elements,elements4,surface_forces,constrained_nodes,volume_force,area_force,u_d,point_forces=point_forces)
u3 = u.reshape((int(len(u)/3),3))
d = ((u3**2).sum(1))**.5
d_max = d.max()
print(d_max)
ax = fea.show3(elements,elements4,triangles_outer,coordinates,u,material,factor=factor) 

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#fea.plot_nodes(coordinates,point_force_nodes,u,ax,factor)
#xyz = fea.compute_deformation(coordinates,u,factor)
#idealab_tools.matplotlib_tools.equal_axes(ax,xyz)
