# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:01:19 2017

@author: danaukes
"""
import pygmsh as pg
import numpy

def convert_2d_to_3d(vertices,z = 0):
    vertices = numpy.array(vertices)
    vertices = numpy.c_[vertices,vertices[:,0:1]*z]
    return vertices

def gen_verts(loop):
    poly4= list(loop.coords)
    poly4 = poly4[:-1]
    poly4 = convert_2d_to_3d(poly4)
    return poly4

def shapely_to_pygmsh(poly,lcar = 1e-1):
    geom = pg.built_in.Geometry()
    
    holes = []
    for interior in poly.interiors:
        hole = geom.add_polygon(gen_verts(interior)[::1], lcar,make_surface=False)
        holes.append(hole.line_loop)
    
    poly = geom.add_polygon(gen_verts(poly.exterior)[::-1],lcar,holes=holes)    
    
    return geom,poly

def build_rgba(A,cmap,normalize=True,invert=False):
    AvC3 = A.copy()
    if invert:
        AvC3 *= -1

    if normalize:
        AvC3 -= AvC3.min()
        AvC3 /= AvC3.max()
    colors = cmap(AvC3)
    colors = colors.squeeze()
    return colors
    