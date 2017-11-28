#!/usr/bin/env python

import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo
import scipy
import cPickle
import matplotlib as mpl
import sivel
import flip

f = open('fix3_x1.pkl','rb')
(fit,_) = pickle.load(f)
ref  = numpy.median(fit['ev_sig'])


mc=[]
f = open('fix3_x1_sim_null0.pkl','rb')
(fit,_) = pickle.load(f)


mc.append(numpy.median(fit['ev_sig']))
for index in xrange(1,100):
    f = open('fix3_x1_sim_null{}.pkl'.format(index),'rb')
    (fit,_) = pickle.load(f)
    mc.append(numpy.median(fit['ev_sig']))
mc= numpy.array(mc)

print ref, mc.max()