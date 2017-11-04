#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel
import matplotlib.pyplot as plt
import corner

# provide stan 3 eigenvectors with respect to try to align delta vector
f = open('temp11.pkl','rb')
(fit,_) = pickle.load(f)

gamma0 = numpy.median(fit['gamma'],axis=0)
gamma1 = numpy.median(fit['rho1'],axis=0)
gamma0_cov = numpy.cov(fit['gamma'],rowvar=False)
gamma1_cov = numpy.cov(fit['rho1'],rowvar=False)


_, gamma0_ev = numpy.linalg.eigh(gamma0_cov)

# print gamma0_ev[:,0]

gamma0_eval = numpy.array(fit['gamma'])
for i in xrange(gamma0_eval.shape[0]):
   for j in xrange(5):
      gamma0_eval[i,j]=numpy.dot(fit['gamma'][i,:], gamma0_ev[:,j])

gamma0_min = numpy.min(gamma0_eval,axis=0)-0.05
gamma0_max = numpy.max(gamma0_eval,axis=0)+0.05

# for i in xrange(5):
#    plt.hist(gamma0_eval[:,i])
#    plt.show()

_, gamma1_ev = numpy.linalg.eigh(gamma1_cov)

gamma1_eval = numpy.array(fit['rho1'])
for i in xrange(gamma1_eval.shape[0]):
   for j in xrange(5):
      gamma1_eval[i,j]=numpy.dot(fit['rho1'][i,:], gamma1_ev[:,j])

gamma1_min = numpy.min(gamma1_eval,axis=0)-0.05
gamma1_max = numpy.max(gamma1_eval,axis=0)+0.05


gamma0median = numpy.median(gamma0_eval,axis=0)
gamma1median = numpy.median(gamma1_eval,axis=0)

fit=None

m = numpy.zeros((2,5))
m[0]=gamma0
m[1]=gamma1
m=m.T

q, r = numpy.linalg.qr(m,mode='complete')
q=q.T

# two color parameter model

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel,sivel_err,x1,x1_err,zcmb,zerr = sivel.sivel(data)

# dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

# dic_meta=cPickle.load(open("META.pkl"))

# sivel=[]
# sivel_err=[]


# for sn in data['snlist']:
#    sn = 'PTF12iiq'
#    if sn in dic_meta.keys() and sn in dic_phreno.keys():
#       meta = dic_meta[sn]
#       vSiII_6355_lbd=0.
#       vSiII_6355_lbd_err=0.
#       counter  = 0
#       for sp in dic_phreno[sn]["spectra"]:
#          if sp in meta['spectra'].keys() and  numpy.abs(meta['spectra'][sp]['salt2.phase']) < 2.5 and numpy.isfinite(dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]):
#             vSiII_6355_lbd += dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
#             vSiII_6355_lbd_err += 1/dic_phreno[sn]['spectra'][sp]["phrenology.vSiII_6355_lbd.err"]**2
#             print dic_phreno[sn]["spectra"][sp]["phrenology.vSiII_6355_lbd"]
#             counter +=1
#       if counter !=0:
#          sivel.append(vSiII_6355_lbd / vSiII_6355_lbd_err)
#          sivel_err.append(1./numpy.sqrt(vSiII_6355_lbd_err)) 
#       else:
#          sivel.append(float('nan'))
#          sivel_err.append(float('nan'))
#    else:
#       sivel.append(float('nan'))
#       sivel_err.append(float('nan'))


# sivel = numpy.array(sivel)
# sivel_err = numpy.array(sivel_err)

use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

snname = numpy.array(data['snlist'])[use]

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn
data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov, \
   'sivel_obs': sivel_renorm, 'sivel_err': sivel_err, 'e1': q[2], 'e2':q[3], 'e3':q[4], 'gamma0in':gamma0,'gamma1in':gamma1,'gamma0in_cov':gamma0_cov,'gamma1in_cov':gamma1_cov,\
   'gamma0_ev':gamma0_ev, 'gamma1_ev':gamma1_ev, 'gamma0_min':gamma0_min, 'gamma0_max': gamma0_max, 'gamma1_min':gamma1_min, 'gamma1_max': gamma1_max }

Delta_simplex = numpy.zeros(nsne-1)
# Delta_simplex = numpy.zeros(nsne)+1./nsne
# k_simplex = numpy.zeros(nsne)
R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()

init = [{'EW' : EW_renorm, \
         'sivel': sivel_renorm,\
         'c_raw' : numpy.zeros(5), \
         'alpha_raw' : numpy.zeros(5), \
         'beta_raw' : numpy.zeros(5), \
         'eta_raw' : numpy.zeros(5), \
         'L_sigma_raw': numpy.zeros(5)+0.03*100, \
         'gamma01': gamma0median[0],\
         'gamma02': gamma0median[1],\
         'gamma03': gamma0median[2],\
         'gamma04': gamma0median[3],\
         'gamma05': gamma0median[4],\
         'gamma11': gamma1median[0],\
         'gamma12': gamma1median[1],\
         'gamma13': gamma1median[2],\
         'gamma14': gamma1median[3],\
         'gamma15': gamma1median[4],\
         'mag_int_raw': mag_renorm, \
         'L_Omega': numpy.identity(5), \
         'Delta_unit':R_simplex, \
         'Delta_scale': 15./4, \
         'k_unit': R_simplex, \
         'k1_unit': R_simplex, \
         'R_unit': numpy.zeros(nsne),\
         # 'rho11': 0./5,\
         # 'rho12': 0./5,\
         # 'rho13': 0./5,\
         'rho1': numpy.zeros(5),\
         } \
        for _ in range(8)]

sm = pystan.StanModel(file='gerard31.stan')
# control = {'stepsize':0.1, 'max_treedepth':20}
control = {'stepsize':1, 'max_treedepth':10}
fit = sm.sampling(data=data, iter=5000, chains=8,control=control,init=init, thin=1)


output = open('temp31.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit
