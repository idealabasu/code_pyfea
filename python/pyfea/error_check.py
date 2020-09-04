# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""

from idealab_tools.data_exchange import dat
import os
import numpy

def compare_output(output,folder):
    for key,value in output.items():
        dat_filename = os.path.join(folder,key+'.dat')
        a = compare_matrices(value,dat_filename)
        if a>0:
            raise(Exception('too many errors'))


def compare_matrices(A,filename,directory=None, format = float,tol = 1e-7):
    directory = directory or ''
    full_path = os.path.normpath(os.path.abspath(os.path.join(directory,filename)))
    B=dat.read(full_path,format = format)
    return num_errors(A,B,tol)

def num_errors(a,b,tol=1e-7):
#    a = numpy.array(a)
#    b = numpy.array(b)
    error = numpy.abs(a-b) > ((numpy.abs(a)).max()*tol)
    num_errors = len(error.nonzero()[0])
    return num_errors

def compare_dict(a,b,tol=1e-7):
    errors = 0
    for key,value in a.items():
        errors+=num_errors(a[key],b[key],tol)
    return errors

def error_check(output,filename,generate_new = False):
        
    import yaml
    import numpy
    
    if generate_new:
        with open(filename,'w') as f:
            for key,value in output.items():
                output[key] = output[key].tolist()
            yaml.dump(output,f)
    else:
    
        with open(filename) as f:
            baseline = yaml.load(f,Loader=yaml.FullLoader)
            for key,value in baseline.items():
                baseline[key] = numpy.array(baseline[key])
        
        if compare_dict(output,baseline)>0:
            raise(Exception('too many errors'))