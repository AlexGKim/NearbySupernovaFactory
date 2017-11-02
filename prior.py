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
from matplotlib import rcParams
rcParams['text.usetex'] = True

dirname = 'output25/'


labels = [r'$\hat{U}$',r'$\hat{B}$',r'$\hat{V}$',r'$\hat{R}$',r'$\hat{I}$']
from matplotlib.ticker import FuncFormatter, MaxNLocator
def format_fn2(tick_val, tick_pos):
    if int(tick_val) in numpy.arange(5):
        return labels[int(tick_val)]
    else:
        return ''

mpl.rcParams['font.size'] = 18

files = ['temp25.pkl','temp27.pkl','temp28.pkl','temp29.pkl','temp30.pkl','temp26.pkl']
labe = [r"$\mbox{{{Cauchy}}}(x_0=0.1,\gamma=0.1)$, $\mbox{{{LKJ}}}(\nu=4)$", \
        r'$\mbox{{{Cauchy}}}(x_0=0.1,\gamma=0.1)$, $\mbox{{{LKJ}}}(\nu=2)$', \
        r'$\mbox{{{Cauchy}}}(x_0=0.1,\gamma=0.1)$, $\mbox{{{LKJ}}}(\nu=8)$', \
        r'$\mbox{{{Cauchy}}}(x_0=0.07,\gamma=0.1)$, $\mbox{{{LKJ}}}(\nu=4)$', \
        r'$\mbox{{{Cauchy}}}(x_0=0.13,\gamma=0.1)$, $\mbox{{{LKJ}}}(\nu=4)$', \
        r'$\mbox{{{Cauchy}}}(x_0=0.1,\gamma=0.1)$, diagonal']
fig = plt.figure()
ax = fig.add_subplot(111)

for ind,(f,lab) in enumerate(zip(files,labe)):
    f = open(f,'rb')
    (fit,_) = pickle.load(f)
    (y, ymin, ymax) = numpy.percentile(fit['rho1']/fit['rho1'][:,0][:,None],(50,50-34,50+34),axis=0)
    ax.errorbar(numpy.arange(5)+0.05*(ind-2),y,yerr=[y-ymin,ymax-y],fmt='o',label=lab)

ax.xaxis.set_major_formatter(FuncFormatter(format_fn2))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.axhline(0,linestyle=':')
ax.set_xlabel(r'Band $X$')
ax.set_xlim((-0.5,4.5))
ax.set_ylabel(r'$\frac{\delta_X}{\delta_{\hat{U}}}$')
# ax.set_ylim((-2,2.))
plt.legend(loc=3,fontsize=14)
pp = PdfPages(dirname+'/delta_prior.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


