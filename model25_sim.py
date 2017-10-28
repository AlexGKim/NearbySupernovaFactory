#!/usr/bin/env python

import pickle
import cPickle
import numpy
import pystan
import sivel

tag="_null"
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
for index in xrange(100):
   pkl_file = open('simdata'+tag+'{}.pkl'.format(index), 'r')
   data = pickle.load(pkl_file)
   pkl_file.close()

   nsne = data['D']
   
   data['e1']=q[2]
   data['e2']=q[3]
   data['e3']=q[4]
   data['gamma0in']=gamma0
   data['gamma1in']=gamma1
   data['gamma0in_cov']=gamma0_cov
   data['gamma1in_cov']=gamma1_cov
   data['gamma0_ev']=gamma0_ev
   data['gamma1_ev']=gamma1_ev
   data['gamma0_min']=gamma0_min
   data['gamma0_max']= gamma0_max
   data['gamma1_min']=gamma1_min
   data['gamma1_max']= gamma1_max

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
            'mag_int_raw': data['mag_obs'], \
            'L_Omega': numpy.identity(5), \
            'Delta_unit':R_simplex, \
            'Delta_scale': 15./4, \
            'k_unit': R_simplex, \
            'k1_unit': R_simplex, \
            'R_unit': numpy.zeros(nsne), \
            # 'rho11': 0./5,\
            # 'rho12': 0./5,\
            # 'rho13': 0./5,\
            'rho1': numpy.zeros(5),\
            } \
           for _ in range(4)]



   sm = pystan.StanModel(file='gerard25.stan')
   control = {'stepsize':1}
   fit = sm.sampling(data=data, iter=5000, chains=4,control=control,init=init, thin=1)


   output = open('temp25_sim'+tag+'{}.pkl'.format(index),'wb')
   pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
   output.close()
   print fit
