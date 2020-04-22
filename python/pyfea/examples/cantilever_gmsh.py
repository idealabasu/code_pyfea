# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import pygmsh as pg
import numpy as np
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
plt.ion()
from mpl_toolkits.mplot3d import Axes3D
import shapely
import shapely.geometry as sg

import idealab_tools.fea_tetra.fea as fea
import idealab_tools.matplotlib_tools


p=  sg.Point(0,0)
poly = p.buffer(1,resolution = 8)
poly = sg.Polygon([(0,0),(1,0),(1,.1),(0,.1)])
poly2= list(poly.exterior.coords)
poly2 = poly2[:-1]
poly2 = numpy.r_[poly2]
poly2 = numpy.c_[poly2,poly2[:,0]*0]
geom = pg.built_in.Geometry()

poly = geom.add_polygon(poly2,lcar=0.05)

axis = [0, 0, .1]
theta = 0

geom.extrude(poly,translation_axis=axis,rotation_axis=axis,point_on_axis=[0, 0, 0], angle=theta)

points, cells, point_data, cell_data, field_data = pg.generate_mesh(geom)
triangles_outer = cells['triangle']

material = fea.Material(100000,.3)
factor = 100

coordinates = points[:]
elements = cells['tetra']

elements,coordinates,mapping = fea.reduce_elements(elements,coordinates)
triangles_outer = numpy.array([[mapping[key] for key in list1] for list1 in triangles_outer],dtype = int)

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

filename = 'cantilever_gmsh.yaml'

output = {}
output['x']=x
output['u']=u

import idealab_tools.fea_tetra.error_check as error_check
error_check.error_check(output,filename,generate_new = False)