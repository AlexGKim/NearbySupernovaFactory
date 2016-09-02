#!/usr/bin/env python

import pickle
import cPickle
import numpy


f = open('temp11.pkl','rb')
(fit, _) = pickle.load(f)
f.close()

for key in fit.keys():
    print key, fit[key].min(), fit[key].max()

print fit['Delta'].flatten().std()

dum = numpy.zeros((5,5))
for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
    dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

dum/= len(fit['L_Omega'])
print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)

dum = numpy.zeros((4,4))
for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
    dum= dum+ numpy.dot(trans,numpy.dot(numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T),trans.T))
dum/= len(fit['L_Omega'])

# color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

pars = ['alpha','alpha','beta','beta','eta','eta','gamma','gamma','rho1','rho1','L_sigma']
pars_n = ['\\alpha_{','{\\alpha ','\\beta_{','{\\beta ','\\eta_{','{\\eta ', '{\\frac{\\gamma}','R_{', '{\\frac{\\delta}','R_{\\delta ','\\sigma_{']
sigfig = [4,1,3,2,4,2,2,2,3]
for p,pn, s in zip(pars,pars_n,sigfig):
    print '${}X}}$'.format(pn)
    for i in xrange(5):
        print '&'
        if pn[0] == 'R':
            dum = numpy.percentile(fit[p][:,i]/(fit[p][:,1]-fit[p][:,2])[:,None],(50,50-34,50+34)) 
        elif pn[0] == '{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34)) 
        else :
            dum = numpy.percentile(fit[p][:,i],(50,50-34,50+34))            
        print  '${1:6.{0}f}^{{{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'


filts = ['U','B','V','R','I']
for i in xrange(5):
    for j in xrange(i+1,5):
        print '{}-{} {} {}'.format(filts[i],filts[j], \
            numpy.std((fit['gamma'][:,i]-fit['gamma'][:,j])[:,None]*fit['k']),numpy.std((fit['rho1'][:,i]-fit['rho1'][:,j])[:,None]*fit['R']))
# print numpy.std(temp[0]),numpy.std(temp[1])

ab = fit['gamma'][:,1][:,None]*fit['R'] + fit['rho1'][:,1][:,None]*fit['R']
av = fit['gamma'][:,2][:,None]*fit['R'] + fit['rho1'][:,2][:,None]*fit['R']
rv = av/(ab-av)

print rv.mean(), rv.std()



mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), \
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']).flatten(),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']).flatten()])

dum = numpy.corrcoef(mega)
print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in dum])


pars = ['gamma','rho1']
pars_n = ['R_{','R_{\\delta ']
sigfig = [3,3]
for p,pn, s in zip(pars,pars_n,sigfig):
    print '${}X}}$'.format(pn)
    for i in xrange(5):
        dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34))            
        print  '${1:6.{0}f}^{{{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'



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

for i in xrange(len(sivel)):
    print '{0} & ${1:5.1f} \pm {2:3.1f}$ & ${3:5.1f} \pm {4:3.1f}$& ${7:5.0f} \pm {8:3.0f}$ & ${5[0]:6.2f} \pm {6[0]:6.2f}$ & ${5[1]:6.2f} \pm {6[1]:6.2f}$& ${5[2]:6.2f} \pm {6[2]:6.2f}$& ${5[3]:6.2f} \pm {6[3]:6.2f}$& ${5[4]:6.2f} \pm {6[4]:6.2f}$ \\\\'.format(i, EW_obs[i,0], numpy.sqrt(EW_cov[i,0,0]),
        EW_obs[i,1], numpy.sqrt(EW_cov[i,1,1]), mag_obs[i,:], numpy.sqrt(numpy.diagonal(mag_cov[i,:,:])),sivel[i],sivel_err[i])

