#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import pickle
import cPickle

f = open('MJC_compile_SNdata.pkl','r')
gal= pickle.load(f)

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
         if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase']) < 2.5  and numpy.isfinite(dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]):
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

i = numpy.intersect1d(names, gal.keys(), assume_unique=True)

inda = numpy.zeros(len(i),dtype='int')

mass=[]
emass = []
for j in xrange(len(i)):
    inda[j] = numpy.where(names == i[j])[0]
    emass.append(gal[i[j]]['eMass'])
    mass.append(gal[i[j]]['Mass'])

mass = numpy.array(mass)
emass= numpy.array(emass)

(x, xmin, xmax) = numpy.percentile(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'],(50,50-34,50+34),axis=0)

plt.errorbar(mass,x[inda],xerr=[emass, emass], yerr=[x[inda]-xmin[inda],xmax[inda]-x[inda]],fmt='o')
plt.ylabel(r'$E_\delta(B-V)$')
plt.xlabel(r'$\log{(M_{host}/M_{\odot})}$')


ux = numpy.array([6.,10])
wm = numpy.where(mass < 10)[0]
print r"${:9.3f} \pm {:9.3f}$".format((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).mean(),(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]),(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 
plt.plot(ux, [numpy.sum(x/dx**2)/dx2,numpy.sum(x/dx**2)/dx2],color='black')
plt.plot(ux, [numpy.sum(x/dx**2)/dx2-1/numpy.sqrt(dx2),numpy.sum(x/dx**2)/dx2-1/numpy.sqrt(dx2)],color='red')
plt.plot(ux, [numpy.sum(x/dx**2)/dx2+1/numpy.sqrt(dx2),numpy.sum(x/dx**2)/dx2+1/numpy.sqrt(dx2)],color='red')


ux = numpy.array([10,13])
wm = numpy.where(mass > 10)[0]
print r"${:9.3f} \pm {:9.3f}$".format((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).mean(),(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).std())
(x, xmin, xmax) = numpy.percentile((((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]),(50,50-34,50+34),axis=0)
dx = (xmax-xmin)/2.
dx2 = numpy.sum(1/dx**2)
print r"${:9.4f} \pm {:9.4f}$".format(numpy.sum(x/dx**2)/dx2, 1/numpy.sqrt(dx2)) 
plt.plot(ux, [numpy.sum(x/dx**2)/dx2,numpy.sum(x/dx**2)/dx2],color='black')
plt.plot(ux, [numpy.sum(x/dx**2)/dx2-1/numpy.sqrt(dx2),numpy.sum(x/dx**2)/dx2-1/numpy.sqrt(dx2)],color='red')
plt.plot(ux, [numpy.sum(x/dx**2)/dx2+1/numpy.sqrt(dx2),numpy.sum(x/dx**2)/dx2+1/numpy.sqrt(dx2)],color='red')

plt.xlim((6,13))

pp = PdfPages("output11/childress.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wm = numpy.where(mass < 10)[0]
low = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()
wm = numpy.where(mass > 10)[0]
hig = (((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None])*fit['R'][:,wm]).flatten()

import scipy.stats
ans= scipy.stats.ks_2samp(hig,low)

plt.hist([low,hig],label=['low mass','high mass'],normed=True,range=[-0.05,0.15])
plt.xlabel(r'$E_\delta(B-V)$')
plt.legend()
pp = PdfPages("output11/childress2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



