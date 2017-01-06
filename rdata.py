#!/usr/bin/env python
import numpy

def rdata():
    f = open('datatable.tex', 'r')
    for i in xrange(3):
        f.readline()

    data = [x.split('&')[0:2] for x in f.readlines()]
    data = data[:-2]
    names = [dum[0].translate(None,"\\").strip() for dum in data]
    names = numpy.array(names)
    temp = [dum[1].replace("+-","-") for dum in data]
    dmu = [dum[0:dum.rfind('\\pm')] for dum in temp]
    ddmu = [float(dum[dum.rfind('{')+1:dum.rfind('}')]) for dum in temp]
    dmu = [float(dum.translate(None,'$')) for dum in dmu]
    dmu = numpy.array(dmu)
    ddmu = numpy.array(ddmu)
    ans = dict(zip(names,zip(dmu,ddmu)))

    return ans