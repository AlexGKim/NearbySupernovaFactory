#!/usr/bin/env python

import numpy

def ellipse_perp(m):
    invm = numpy.linalg.inv(m)
    y = numpy.dot([0.,1], numpy.dot(invm,[0,1]))
    y=numpy.sqrt(1/y)
    x = numpy.dot([1.,0], numpy.dot(invm,[1.,0]))
    x=numpy.sqrt(1/x)
    return x,y