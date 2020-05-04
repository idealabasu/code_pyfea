# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:49:07 2017

@author: danaukes
"""

import PyQt5.QtGui as qg
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc

import os
import sys
import pyqtgraph.opengl as pgo

import pyfea_examples.tube.tube_merged as tg
import logging
import traceback

base = os.path.abspath(os.path.normpath(os.path.expanduser('~')))
filename = os.path.join(base,'error_log.txt')

logger = logging.Logger('design_tool',level=logging.DEBUG)
handler = logging.FileHandler(filename=filename,mode='w')
logger.addHandler(handler)  
excepthook_internal = sys.excepthook

def excepthook(exctype,value,tb):
    if exctype is not SystemExit:
        message = '''{}: {}'''.format(str(exctype),str(value))
        print(message)

        tbmessage = traceback.format_tb(tb)
        tbmessage = '  '.join(tbmessage)

        logger.error(message)
        logger.debug('\n'+tbmessage)
        
        excepthook_internal(exctype,value,tb)

class W2(pgo.GLViewWidget):
    nominal_width = 1024
    nominal_height = 768
    def sizeHint(self):
        return qc.QSize(self.nominal_width, self.nominal_height)    



class Widget(qw.QWidget):
    nominal_width = 1024
    nominal_height = 768
    def __init__(self):
        super(Widget,self).__init__()
        
        layout1 = qw.QHBoxLayout()
        self.text = ''
        
        self.w = W2() 
        
        self.w.opts['center'] = qg.QVector3D(0,0,.05)
        self.w.opts['distance'] = .1
        self.w.opts['azimuth'] = -45
        self.w.opts['elevation'] = 30
#        self.w.resize(*size)

        self.button_gen = qw.QPushButton('Generate')
        self.button_comp = qw.QPushButton('Compute')
        self.field_d1_outer = qw.QLineEdit('.4')
        self.field_d2_outer = qw.QLineEdit('.2')
        self.field_d1_inner = qw.QLineEdit('.2')
        self.field_d2_inner = qw.QLineEdit('.1')
        self.characteristic_length_max = qw.QLineEdit('.2')
        self.characteristic_length_min = qw.QLineEdit('.1')
        self.field_length = qw.QLineEdit('1')
        self.field_youngs = qw.QLineEdit('1e6')
        self.field_poisson = qw.QLineEdit('.3')
        self.field_density = qw.QLineEdit('1000')
        self.field_max_stress = qw.QLineEdit()
#        w.addItem(mi)

        self.field_output = qw.QTextEdit()
        layout2 = qw.QVBoxLayout()
        
        layout1_1 = qw.QHBoxLayout()
        layout1_1.addWidget(self.button_gen,1)
        layout1_1.addWidget(self.button_comp)
        layout1_1.addStretch()

        fields = self.field_d1_outer,self.field_d2_outer,self.field_d1_inner,self.field_d2_inner,self.field_length,self.characteristic_length_max,self.characteristic_length_min,self.field_youngs,self.field_poisson,self.field_density
        labels = 'd1_outer','d2_outer','d1_inner','d2_inner','length','characteristic_length_max','characteristic_length_min','youngs', 'poisson','density'
        
        for label,field in zip(labels,fields):
            layout = qw.QVBoxLayout()
            layout.addWidget(qw.QLabel(label))
            layout.addWidget(field)
            layout2.addLayout(layout)
        
        
        layout2.addStretch()
        layout2.addLayout(layout1_1)
        layout2.addWidget(self.field_output)
        
        layout = qw.QVBoxLayout()
        layout.addWidget(qw.QLabel('max_stress'))
        layout.addWidget(self.field_max_stress)
        layout2.addLayout(layout)
        
#        layout1.addStretch()
        layout1.addWidget(self.w)
        layout1.addLayout(layout2)
        self.setLayout(layout1)
        
        self.button_gen.clicked.connect(self.generate_geometry_outer)
        self.button_comp.clicked.connect(self.compute_outer)
        self.field_max_stress.setEnabled(False)


    def write(self, text):
        self.text+=text
        self.field_output.setText(self.text)

    def sizeHint(self):
        buffer_x = 14
        buffer_y = 36
        return qc.QSize(self.nominal_width - buffer_x, self.nominal_height - buffer_y)    
    
    def generate_geometry_outer(self):
        d1_outer = float(self.field_d1_outer.text())
        d2_outer = float(self.field_d2_outer.text())
        d1_inner = float(self.field_d1_inner.text())
        d2_inner = float(self.field_d2_inner.text())
        length = float(self.field_length.text())
        characteristic_length_max = float(self.characteristic_length_max.text())
        characteristic_length_min = float(self.characteristic_length_min.text())
        coordinates, triangles_outer,elements  = tg.generate_geometry(length,d1_outer,d2_outer,d1_inner,d2_inner,char_len_max=characteristic_length_max,char_len_min=characteristic_length_min)
        self.coordinates = coordinates
        self.triangles_outer = triangles_outer
        self.elements = elements  

    def compute_outer(self):
        self.field_output.clear()
        for item in self.w.items:
            self.w.removeItem(item)
        poisson = float(self.field_poisson.text())
        youngs = float(self.field_youngs.text())
        density = float(self.field_density.text())
        
        mi,c_max = tg.compute(self.coordinates, self.triangles_outer,self.elements,youngs,poisson,density)

        self.w.addItem(mi)
        print('max stress: ',c_max)
        self.field_max_stress.setText(str(c_max))


if __name__=='__main__':        
    app = qg.QApplication(sys.argv)
    sys.excepthook = excepthook          


    w = Widget()    

    sys.stdout = w
    sys.stderr = w

    w.show()
    sys.exit(app.exec_())
