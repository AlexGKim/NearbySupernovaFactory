#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle

f = open('lsSFR.txt', 'r')
for i in xrange(1):
    f.readline()

data = [x.split() for x in f.readlines()]
f.close()

rnames=[]
lsSFR = []
mlsSFR = []
plsSFR = []


for d in data:
    rnames.append(d[0].translate(None, "'[(,])'"))
    lsSFR.append(float(d[1].translate(None, '[(,])')))
    mlsSFR.append(float(d[3].translate(None, '[(,])')))
    plsSFR.append(float(d[2].translate(None, '[(,])')))


rnames = numpy.array(rnames)
lsSFR = numpy.array(lsSFR)
mlsSFR=numpy.array(mlsSFR)
plsSFR = numpy.array(plsSFR)


f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

dic_meta=cPickle.load(open("META.pkl"))

sivel=[]
sivel_err=[]
for sn in data['snlist']:
   if sn in dic_meta.keys() and sn in dic_phreno.keys():
      meta = dic_meta[sn]
      vSiII_6355_lbd=0.
      vSiII_6355_lbd_err=0.
      counter  = 0
      for sp in dic_phreno[sn]["spectra"]:
         if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase'] < 3):
            vSiII_6355_lbd += dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
            vSiII_6355_lbd_err += 1/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
            counter +=1
      if counter !=0:
         sivel.append(vSiII_6355_lbd / vSiII_6355_lbd_err)
         sivel_err.append(1./numpy.sqrt(vSiII_6355_lbd_err))
      else:
         sivel.append(float('nan'))
         sivel_err.append(float('nan'))
   else:
      sivel.append(float('nan'))
      sivel_err.append(float('nan'))

sivel = numpy.array(sivel)
sivel_err = numpy.array(sivel_err)

use = numpy.isfinite(sivel)

names= numpy.array(data['snlist'])[use]

i = numpy.intersect1d(names, rnames, assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')
indr = numpy.zeros(len(i),dtype='int')

for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    indr[j] = numpy.where(rnames == i[j])[0]

(x, xmin, xmax) = numpy.percentile(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'],(50,50-34,50+34),axis=0)
plt.errorbar(lsSFR[indr],x[inda],xerr=[mlsSFR[indr], plsSFR[indr]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.ylabel(r'$E_\delta(B-V)$')
plt.xlabel(r'lsSFR')
pp = PdfPages("output11/lsSFR.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.show()