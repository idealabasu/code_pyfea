# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

import meshio
points,cells,a,b,c= meshio.read('Part2.off')
meshio.write('Part2.msh', points, cells)

