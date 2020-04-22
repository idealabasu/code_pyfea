# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:18:44 2017

@author: danaukes
"""

import PyQt5.QtGui as qg
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc

import sys
import pyqtgraph.opengl as pgo

#import pygmsh as pg
import pygmsh
import numpy as np
import numpy
#        import matplotlib.cm as cm
#        import matplotlib.pyplot as plt
#        plt.ion()
#        from mpl_toolkits.mplot3d import Axes3D
#import shapely
#import shapely.geometry as sg

import pyfea.fea as fea
#        import idealab_tools.matplotlib_tools
#import idealab_tools.pygmsh_fix

import pyfea.fea as fea
import pygmsh
import idealab_tools.plot_tris

def generate_geometry(l,d1_outer,d2_outer,d1_inner,d2_inner,char_len_min = .1,char_len_max = .2):
    geom = pygmsh.opencascade.Geometry(characteristic_length_min=char_len_min,characteristic_length_max=char_len_max)
    
    c1 = geom.add_cone([0,0,0], [0,0,l],d1_outer,d1_inner)
    c2 = geom.add_cone([0,0,0], [0,0,l],d2_outer,d2_inner)
    
    c3 = geom.boolean_difference([c1], [c2])
    
#    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
    mo =pygmsh.generate_mesh(geom)
    triangles_outer = mo.cells['triangle']
    #face_colors = [(1,0,0,1)]*len(triangles_outer)
    #face_colors = numpy.array(face_colors)
    #idealab_tools.plot_tris.plot_tris(points,triangles_outer,face_colors = face_colors,drawEdges=True, edgeColor = (1,1,1,1))
#    youngs,poisson=1e5,1e-2
    

    
    coordinates = mo.points[:]
    elements = mo.cells['tetra']
    
    used_elements = fea.find_used_elements(elements,triangles_outer)
    coordinates,mapping = fea.coord_reduce(coordinates,used_elements)
    triangles_outer = fea.element_reduce(triangles_outer,mapping)
    elements= fea.element_reduce(elements,mapping)
#    del points
    
    a=fea.analyze(coordinates,elements)
    print(a)
    elements[a] = elements[a][:,(0,2,1,3)]
    a=fea.analyze(coordinates,elements)
    print(a)
    
    T = coordinates[elements[:,1:]]-coordinates[elements[:,0:1]]
    dt = numpy.array([numpy.linalg.det(item) for item in T])
    elements = elements[dt!=0]
    return coordinates, triangles_outer , elements,  

def compute(coordinates, triangles_outer, elements,youngs,poisson,density,lcar = 1e-1):

    material = fea.Material(youngs,poisson)
    factor = 1
#    f = 1e-2    
    
    xx = coordinates[:,0]
    yy = coordinates[:,1]
    zz = coordinates[:,2]
    
    z_max = coordinates.max(0)[2]
    z_min = coordinates.min(0)[2]
    x_max = coordinates.max(0)[0]
    x_min = coordinates.min(0)[0]
    
    
    
    ii_bottom = ((coordinates[triangles_outer,2]==z_min).sum(1)==3)
    ii_top = ((coordinates[triangles_outer,2]==z_max).sum(1)==3)
    dirichlet_bottom = triangles_outer[ii_bottom]
    dirichlet_top = triangles_outer[ii_top]
    #dirichlet = numpy.r_[dirichlet_bottom,dirichlet_top]
    dirichlet = dirichlet_bottom
    
    neumann = numpy.zeros((0,3),dtype = int)
    
    dirichlet_nodes = numpy.unique(dirichlet)
    neumann_nodes = numpy.unique(neumann)
    
    #fea.plot_triangles(coordinates,triangles_outer)
#    density = 1000
    
    def volume_force(x):
        volforce = numpy.zeros((x.shape[0],3))
        volforce[:,1] = 9.81*density
        return volforce    
    
    def u_d(x):
        mm = x.shape[0]
        M = numpy.zeros((3*mm,3))
        W = numpy.zeros((3*mm,1))
        M[0::3,0] = 1
        M[1::3,1] = 1
        M[2::3,2] = 1
    
    #    ii = (x[:,2]==z_max).nonzero()[0]
    #    W[ii*3+2] = f * x[ii,0:1]+f
        return W,M
    
    u_max = []
    C_max = []
    x,u = fea.compute(material,coordinates,elements,[],neumann,dirichlet_nodes,volume_force,fea.surface_force_empty,u_d)
    C = fea.max_stress(elements,coordinates, u,material)
    u_max.append([u[0::3].max(),u[2::3].max(),u[2::3].max()])
    C_max.append(C.max())
    
    filename = 'tube.yaml'

    output = {}
    output['x']=x
    output['u']=u

#    import pyfea.error_check as error_check
#    error_check.error_check(output,filename,generate_new = False)

    xyz,triangles_outer,c_face,c_vertex = fea.show23_int(elements,triangles_outer,coordinates,u,material,factor)

    md = pgo.MeshData(vertexes = xyz,faces = triangles_outer,vertexColors = c_vertex)
    mi = pgo.GLMeshItem(meshdata = md,shader='balloon',drawEdges=False,edgeColor = [1,1,1,1],smooth=False,computeNormals = False,glOptions='opaque')
    return mi,C.max()
#        ax = fea.show3(elements,triangles_outer,coordinates,u,material,factor = factor)

if __name__=='__main__':
    coordinates, triangles_outer,elements  = generate_geometry(1,.4,.2,.2,.1)
    mi,C_max = compute(coordinates, triangles_outer,elements,1e5,1e-2,1000)
    import idealab_tools.plot_tris as pt
    pt.plot_mi(mi)