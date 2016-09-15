#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle

f = open('forAlex_tmp14sept2016_localhost_parameters.dat', 'r')
for i in xrange(1):
    f.readline()

data = [x.split() for x in f.readlines()]
f.close()

rnames=[]
lsSFR = []
mlsSFR = []
plsSFR = []
mass = []
pmass = []
mmass = []
sfr = []
psfr=[]
msfr=[]
pdelay =[]

for d in data:
    rnames.append(d[0])
    lsSFR.append(float(d[1]))
    mlsSFR.append(float(d[3]))
    plsSFR.append(float(d[2]))
    mass.append(float(d[4]))
    mmass.append(float(d[6]))
    pmass.append(float(d[5]))
    sfr.append(float(d[7]))
    msfr.append(float(d[9]))
    psfr.append(float(d[8]))
    pdelay.append(float(d[10]))

rnames = numpy.array(rnames)
lsSFR = numpy.array(lsSFR)
mlsSFR=numpy.array(mlsSFR)
plsSFR = numpy.array(plsSFR)
mass = numpy.array(mass)
mmass=numpy.array(mmass)
pmass = numpy.array(pmass)
sfr = numpy.array(sfr)
msfr=numpy.array(msfr)
psfr = numpy.array(psfr)
pdelay = numpy.array(pdelay)

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
fig, axes = plt.subplots(nrows=4)
axes[0].errorbar(lsSFR[indr],x[inda],xerr=[mlsSFR[indr], plsSFR[indr]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
axes[0].set_ylabel(r'$E_\delta(B-V)$')
axes[0].set_xlabel(r'lsSFR')
axes[1].errorbar(mass[indr],x[inda],xerr=[mmass[indr], pmass[indr]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
axes[1].set_ylabel(r'$E_\delta(B-V)$')
axes[1].set_xlabel(r'mass')
axes[2].errorbar(sfr[indr],x[inda],xerr=[msfr[indr], psfr[indr]], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
axes[2].set_ylabel(r'$E_\delta(B-V)$')
axes[2].set_xlabel(r'sfr')
axes[3].errorbar(pdelay[indr],x[inda], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
axes[3].set_ylabel(r'$E_\delta(B-V)$')
axes[3].set_xlabel(r'pdelay')
fig.subplots_adjust(hspace=.3)
fig.set_size_inches(8,11)
pp = PdfPages("output11/rigault.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.show()