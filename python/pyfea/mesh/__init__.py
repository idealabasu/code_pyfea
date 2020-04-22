# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""
import pygmsh as pg
import numpy
import yaml
import os
from idealab_tools.data_exchange.generic_data import GenericData

def generate_plate():
    geom = pg.opencascade.Geometry(characteristic_length_max=.0025, characteristic_length_min=.001)
    geom.add_box((0,0,0),(.10,.10,.005))
    mo = pg.generate_mesh(geom)
    return mo



def mesh_to_dict(mo):
    mo2 = {}
    mo2['cell_data'] = mo.cell_data
    mo2['cells'] = {}
    for key,value in mo.cells.items():
        mo2['cells'][key] = mo.cells[key].tolist(),mo.cells[key].dtype
    mo2['field_data'] = mo.field_data
    mo2['gmsh_periodic'] = mo.gmsh_periodic
    mo2['info'] = mo.info
    mo2['node_sets'] = mo.node_sets
    mo2['point_data'] = mo.point_data
    mo2['points'] = mo.points.tolist(),mo.points.dtype
    return mo2    

def save_mesh(mo,meshname):
    #mo2['prune']
    mo2 = mesh_to_dict(mo)
    
    with open(meshname,'w') as f:
        yaml.dump(mo2,f)
    
def load_mesh(meshname):
    with open(meshname) as f:
        mo2 = yaml.load(f)
    cells = {}
    for key,value in mo2['cells'].items():
        cells[key] = numpy.array(*value)
    mo2['cells'] = cells
    mo2['points'] = numpy.array(*(mo2['points']))
    mo = GenericData(**mo2)
    return mo


def get_mesh(meshname,generator,force_regen = False):
    if not force_regen:
        if os.path.exists(meshname):
            try:
                mo = load_mesh(meshname)
            except Exception:
                mo=generator()
                save_mesh(mo,meshname)
        else:
            mo=generator()
            save_mesh(mo,meshname)
    else:
        mo=generator()
        save_mesh(mo,meshname)
    return mo