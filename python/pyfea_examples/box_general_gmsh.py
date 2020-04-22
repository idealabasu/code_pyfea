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

import pyfea.fea as fea
import idealab_tools.matplotlib_tools

def select_tri_from_nodes(points,tris,points_filter):
    return points_filter[tris].sum(1)==3

def select_quad_from_nodes(points,quads,points_filter):
    return points_filter[quads].sum(1)==4

def create_t(points,quads):
    quad_points = points[quads,:]
    dqp = quad_points[1:] - quad_points[0]

geom = pg.opencascade.Geometry(characteristic_length_max=.1, characteristic_length_min=.1)
geom.add_box((0,0,0),(1,1,1))

points, cells, point_data, cell_data, field_data = pg.generate_mesh(geom)
triangles_outer = cells['triangle']

material = fea.Material(100000,.3)
factor = 100

coordinates = points[:]
elements = cells['tetra']

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

z_max = coordinates.max(0)[2]
z_min = coordinates.min(0)[2]
x_max = coordinates.max(0)[0]
x_min = coordinates.min(0)[0]



ii_bottom = ((coordinates[triangles_outer,0]==x_min).sum(1)==3)
ii_top = ((coordinates[triangles_outer,0]==x_max).sum(1)==3)
#ii_neumann = (ii_bottom+ii_top)==0
dirichlet_bottom = triangles_outer[ii_bottom]
dirichlet_top = triangles_outer[ii_top]
dirichlet = numpy.r_[dirichlet_bottom,dirichlet_top]

neumann = numpy.zeros((0,3),dtype = int)

dirichlet_nodes = numpy.unique(dirichlet)
neumann_nodes = numpy.unique(neumann)

#heat_source_nodes = 
#fea.plot_triangles(coordinates,triangles_outer)

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

x,u = fea.compute(material,coordinates,elements,neumann,dirichlet_nodes,fea.volume_force_empty,fea.surface_force_empty,u_d)
ax = fea.show(elements,triangles_outer,coordinates,u,material,factor=factor) 
fea.plot_nodes(coordinates,dirichlet_nodes,u,ax,factor)