#!/usr/bin/env python
import pickle
import numpy
import pystan
import sncosmo
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

synname=['U','B','V','R','I']


snmod='hsiao'

synlam = numpy.array([[3300.00, 3978.02]
    ,[3978.02,4795.35]
    ,[4795.35,5780.60]
    ,[5780.60,6968.29]
    ,[6968.29,8400.00]])

synname=['U','B','V','R','I']

synbands=[]

for name, lams in zip(synname,synlam):
    synbands.append(sncosmo.Bandpass(lams, [1.,1.], name='tophat'+name))

model_nodust = sncosmo.Model(source=snmod)
flux_nodust = model_nodust.bandflux(synbands,0.)

a = [1.]
rv = numpy.exp(numpy.arange(numpy.log(.9), numpy.log(8)+0.001,numpy.log(8/.9)/50))

avs=[]
ebvs=[]
rvs=[]
AX = []

for r in rv:
    dust = sncosmo.F99Dust(r_v=r)
    dust.set(ebv=a/r)
    model = sncosmo.Model(source=snmod, effects=[dust], effect_names=['host'], effect_frames=['rest'])
    dust2 = sncosmo.CCM89Dust()
    dust2.set(ebv=a/r,r_v=r)
    model2 = sncosmo.Model(source=snmod, effects=[dust2], effect_names=['host'], effect_frames=['rest'])
    AX.append(-2.5*numpy.log10(model.bandflux(synbands,0.)/flux_nodust)+2.5*numpy.log10(model2.bandflux(synbands,0.)/flux_nodust))
    avs.append(a)
    ebvs.append(a/r)
    rvs.append(r)

avs = numpy.array(avs)
ebvs = numpy.array(ebvs)
AX = numpy.array(AX)
rvs=numpy.array(rvs)


for i in xrange(5):
    plt.plot(rvs,AX[:,i],label=synname[i])
plt.xlabel(r'$R$')
plt.ylabel(r'$A^C_X-A^F_X$ @ $A_V=1$')
plt.legend(loc=4)
pp = PdfPages('output18/CCMF99diff.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()