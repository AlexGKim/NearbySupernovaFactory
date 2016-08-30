#!/usr/bin/env python

import pickle
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

pars = ['alpha','beta','eta','gamma','rho1','L_sigma']
pars_n = ['\\alpha_{','\\beta_{','\\eta_{','R_{','R_{\\delta ','\\sigma_{']
for p,pn in zip(pars,pars_n):
    print '${}X}}$'.format(pn)
    for i in xrange(5):
        print '&'

        if pn == 'R_{' or pn == 'R_{\\delta ':
            dum = numpy.percentile(fit[p][:,i]/(fit[p][:,1]-fit[p][:,2])[:,None],(50,50-34,50+34)) 
        else :
            dum = numpy.percentile(fit[p][:,i],(50,50-34,50+34))            
        print  '${:6.4f}^{{{:6.4f}}}_{{{:6.4f}}}$'.format(dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'


temp = numpy.array([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['R'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']])
print numpy.std(temp[0]),numpy.std(temp[1])

ab = fit['gamma'][:,1][:,None]*fit['R'] + fit['rho1'][:,1][:,None]*fit['R']
av = fit['gamma'][:,2][:,None]*fit['R'] + fit['rho1'][:,2][:,None]*fit['R']
rv = av/(ab-av)
print rv.mean(), rv.std()