# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 18:32:31 2013

@author: bjchen
"""

import numpy
import inspect

def savedata(filename, *args):
    stack = inspect.stack()
    try:
        locals_ = stack[1][0].f_locals
        varstr = ','.join(['%s=locals_[\'%s\']'%(v,v) for v in args])        
        eval('numpy.savez_compressed(filename, %s)'%varstr) 
    finally:
        del stack
    
def loaddata(filename):
    stack = inspect.stack()
    tmp = numpy.load(filename)
    try:
        locals_ = stack[1][0].f_locals        
        for k in tmp.keys():
            locals_[k] = tmp[k]
    except:
        del stack
        tmp.close()
    tmp.close()
    
def ismember(a, b):
    if type(b) is not numpy.ndarray: 
        b = numpy.array(b)
    tf = numpy.in1d(a,b) # for newer versions of numpy
    index = numpy.array([(numpy.where(b == i))[0][-1] if t else 0 for i,t in zip(a,tf)])
    return tf, index
    
def unionrows(a, b, ai=False, bi=False):
    if type(a) is not numpy.ndarray:
        a = numpy.array(a)
    if type(b) is not numpy.ndarray:
        b = numpy.array(b)
    nr, nc = a.shape
    datatype = {'names':['f{}'.format(i) for i in range(nc)], 'formats':nc*[a.dtype] }
    tmpu = numpy.union1d(a.view(datatype), b.view(datatype))
    u = tmpu.view(a.dtype).reshape(-1, nc)
    if not ai and not bi:
        return u
        
    if ai:
        tf, aidx = ismember(a.view(datatype), tmpu)
    if bi:
        tf, bidx = ismember(b.view(datatype), tmpu)
    
    if ai and bi:
        return u, aidx, bidx
    elif ai and not bi:
        return u, aidx
    elif (not ai) and bi:
        return u, bidx
        