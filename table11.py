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

pars = ['alpha','alpha','beta','beta','eta','eta','gamma','rho1','L_sigma']
pars_n = ['\\alpha_{','{\\alpha ','\\beta_{','{\\beta ','\\eta_{','{\\eta ','R_{','R_{\\delta ','\\sigma_{']
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
