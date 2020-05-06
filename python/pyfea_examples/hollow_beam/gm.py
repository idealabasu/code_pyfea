# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import os
import subprocess
import meshio
import numpy
import pygmsh as pg
import numpy as np
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yaml
from idealab_tools.data_exchange.generic_data import GenericData

import pyfea.fea as fea
#
##def makeMesh(filename, delete_files=True):
#template = '''
#Geometry.LineNumbers = 0;
#Geometry.SurfaceNumbers = 0;
#Merge "{0}";
#Surface Loop(1) = {{1}};
#//+
#Volume(1) = {{1}};
#'''
#
#geo_string = template.format('beam.STL')
#with open('output.geo', 'w') as f:
#    f.writelines(geo_string)
#
#command_string = 'gmsh output.geo -3 -format msh'
#p = subprocess.Popen(command_string, shell=True)
#p.wait()
#mesh_file = 'output.msh'
#data = meshio.read(mesh_file)
##if delete_files:
##    os.remove('output.msh')
##    os.remove('output.geo')
##return data



import meshio
mo = meshio.read('output.msh')

points = mo.points
#tet = mo.cells['tetra']
#tris = mo.cells['triangle']

triangles_outer = mo.cells['triangle']
coordinates = mo.points[:]
elements = mo.cells['tetra']


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

def volume_force(x):
    density = 1000
    volforce = numpy.zeros((x.shape[0],3))
    volforce[:,1] = -9.81*density
#    volforce[:,2] = -9.81
    return volforce    

def area_force(x,n):
    f = numpy.zeros((1,3))
#    f = numpy.array([[0,0,-p]])
    return f


point_forces = []

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

material = fea.Material(1e8,.3)
factor = 1

x,u = fea.compute(material,coordinates,elements,elements4,surface_forces,constrained_nodes,volume_force,area_force,u_d,point_forces=point_forces)
u3 = u.reshape((int(len(u)/3),3))
d = ((u3**2).sum(1))**.5
d_max = d.max()
print(d_max)
fea.show3(elements,elements4,triangles_outer,coordinates,u,material,factor=factor,draw_edges = True) 
