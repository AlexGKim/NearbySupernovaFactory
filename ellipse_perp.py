#!/usr/bin/env python

import numpy

def anglestuff(l1,l2,a1,a2):
    if a2==0:
        cosphi=0.
        sinphi=1.
    else:
        temp = numpy.arctan(-l1*a1/l2/a2)
        cosphi = numpy.cos(temp)
        sinphi = numpy.sin(temp)  
    return cosphi,sinphi

def ellipse_perp(m):
    # 1, 2, eigenvector; a,b == x,y 
    (l1, l2),((a1, b1),(a2,b2)) = numpy.linalg.eig(m)

    cosphi,sinphi = anglestuff(l1,l2,a1,a2)
    y = l1*b1*cosphi + l2*b2*sinphi

    cosphi,sinphi = anglestuff(l1,l2,b1,b2)
    x = l1*a1*cosphi + l2*a2*sinphi
    return numpy.abs(x),numpy.abs(y)