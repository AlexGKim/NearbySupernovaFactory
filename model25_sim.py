#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel
f = open('temp11.pkl','rb')
(fit,_) = pickle.load(f)

gamma0 = numpy.median(fit['gamma'],axis=0)
gamma1 = numpy.median(fit['rho1'],axis=0)
fit=None

m = numpy.zeros((2,5))
m[0]=gamma0
m[1]=gamma1
m=m.T

q, r = numpy.linalg.qr(m,mode='complete')
q=q.T


# two color parameter model
for index in xrange(1):
   pkl_file = open('simdata_null{}.pkl'.format(index), 'r')
   data = pickle.load(pkl_file)
   pkl_file.close()

   nsne = data['D']
   
   data['e1']=q[2]
   data['e2']=q[3]
   data['e3']=q[4]

   Delta_simplex = numpy.zeros(nsne-1)
   # Delta_simplex = numpy.zeros(nsne)+1./nsne
   # k_simplex = numpy.zeros(nsne)
   R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
   R_simplex = R_simplex/R_simplex.sum()

   init = [{'EW' : data['EW_obs'], \
            'sivel': data['sivel_obs'],\
            'c_raw' : numpy.zeros(5), \
            'alpha_raw' : numpy.zeros(5), \
            'beta_raw' : numpy.zeros(5), \
            'eta_raw' : numpy.zeros(5), \
            'L_sigma_raw': numpy.zeros(5)+0.03*100, \
            'gamma01': 61./5,\
            'gamma02': 50./5,\
            'gamma03': 40./5,\
            'gamma04': 30./5,\
            'gamma05': 20./5,\
            'gamma11': -12./5,\
            'gamma12': -14./5,\
            'gamma13': -16./5,\
            'gamma14': -14./5,\
            'gamma15': -14/5,\
            'mag_int_raw': data['mag_obs'], \
            'L_Omega': numpy.identity(5), \
            'Delta_unit':R_simplex, \
            'Delta_scale': 15./4, \
            'k_unit': R_simplex, \
            'k1_unit': R_simplex, \
            'R_unit': R_simplex, \
            'rho11': 0./5,\
            'rho12': 0./5,\
            'rho13': 0./5,\
            } \
           for _ in range(4)]



   sm = pystan.StanModel(file='gerard25.stan')
   control = {'stepsize':1}
   fit = sm.sampling(data=data, iter=5000, chains=4,control=control,init=init, thin=1)


   output = open('temp25_sim_null{}.pkl'.format(index),'wb')
   pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
   output.close()
   print fit