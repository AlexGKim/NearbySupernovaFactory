#!/usr/bin/env python

import pickle
import cPickle
import numpy
import sivel
import sncosmo

f = open('temp23.pkl','rb')
(fit, _) = pickle.load(f)
f.close()


for key in fit.keys():
    print key, fit[key].min(), fit[key].max()

print numpy.std(((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + \
  (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1']).flatten())

print numpy.std(((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']).flatten())



# # print 'fit data to Fitzpatrick'
# gamma =  fit['gamma']
# gammacov = numpy.cov(gamma,rowvar=False)
# gamma = numpy.mean(gamma,axis=0)
# gammainvcov = numpy.linalg.inv(gammacov)

# def lnprob(p, x, y ,yerr):
#   # p is ebv and r_v
#   f99 = sncosmo.F99Dust(r_v =p[1])
#   f99.set(ebv=p[0])
#   A_ = f99.propagate(x,1.)
#   A_ = -2.5*numpy.log10(A_)
#   dum = y-A_
#   ans = -0.5*numpy.dot(dum, numpy.dot(yerr, dum))
#   return ans

# efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])
# ndim, nwalkers = 2, 100
# pos = [numpy.array([12., 2.5]) + 1e-4*numpy.random.randn(ndim) for i in range(nwalkers)]
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(efflam, gamma, gammainvcov))
# sampler.run_mcmc(pos,2000)
# samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
# (y,ymin,ymax)  = numpy.percentile(samples[:,1],(50,50-34,50+34)) 

# print 'Fitzpatrick fit'
# print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)
# (y,ymin,ymax)  = numpy.percentile(samples[:,0],(50,50-34,50+34)) 
# print '{:6.2f} {:6.2f} {:6.2f}'.format(y, y-ymin, ymax-y)

print "probablility of no intrinsic"
# print numpy.min(fit['rho1'],axis=0)
# print numpy.max(fit['rho1'],axis=0)


mod = numpy.sqrt(numpy.sum(fit['rho1']**2,axis=1))
bin=0.5
ran=((-.98,5.3+bin),(-4.5,1.7+bin),(-4.3,1.3+bin),(-2,2.3+bin),(0,5.8+bin))
edges=[];
for i in ran:
  edges.append(numpy.arange(i[0],i[1]+0.00001,bin))

H,edges= numpy.histogramdd(fit['rho1'],bins=edges)
print H.max()

print H.shape, H.sum()
indeces=numpy.zeros(5,dtype='int')
for ind in xrange(5):
  indeces[ind] = numpy.digitize(0,edges[ind])

indeces=indeces-1
print indeces
print H.sum()
print H[indeces[0],indeces[1],indeces[2],indeces[3],indeces[4]]
wbig = H > H[indeces[0],indeces[1],indeces[2],indeces[3],indeces[4]]
print 'confidence non-zero ',H[wbig].sum() / H.sum()
print 'confidence non-zero ',1./(H >0).sum(), (H >0).sum()

print "standard deviation of delta"
print fit['Delta'].flatten().std()

dum = numpy.zeros((5,5))
for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
    dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

dum/= len(fit['L_Omega'])
print "average L_Omega"
print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)

dum = numpy.zeros((4,4))
for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
    dum= dum+ numpy.dot(trans,numpy.dot(numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T),trans.T))
dum/= len(fit['L_Omega'])

# color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

print "color covariance"
dumsig = numpy.sqrt(numpy.diag(dum))
print [" , ".join(map('{0:.3f}'.format, dumsig))]
dumcor =  dum/ numpy.dot(dumsig[:,None],dumsig[None,:])
print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in dumcor])

print "standard deviations of E(B-V)"
print numpy.std(fit['k']*((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None])
print numpy.std(fit['R']*((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None])

pars = ['gamma', 'gamma1','rho1']
pars_n = ['\\gamma^0_X',  '\\gamma^1_X', '\\delta_X']
sigfig = [3,3,3]
for p,pn, s in zip(pars,pars_n,sigfig):
    print '${}$'.format(pn)
    for i in xrange(5):
        print '&'
        if pn[0] == 'R':
            dum = numpy.percentile(fit[p][:,i]/(fit[p][:,1]-fit[p][:,2])[:,None],(50,50-34,50+34)) 
        elif pn[0:2] == '{{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,4]-1,(50,50-34,50+34)) 
        elif pn[0] == '{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34)) 
        else :
            dum = numpy.percentile(fit[p][:,i],(50,50-34,50+34))            
        print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'


mega = numpy.array([fit['k'].flatten(),fit['k1'].flatten(),fit['R'].flatten()])
dum = numpy.corrcoef(mega)
print "k-D correlation coefficients"
print " \\\\\n".join([" & ".join(map('{0:.2f}'.format, line)) for line in dum])



pars = ['alpha','alpha','beta','beta','eta','eta','gamma','gamma','gamma1','gamma1','rho1','rho1','L_sigma']
pars_n = ['\\alpha_X','{\\alpha_X/\\alpha_V-1}','\\beta_X','{\\beta_X/\\beta_V-1}',\
  '\\eta_X','{\\eta_X/\\eta_V-1}', '\\gamma^0_X', '{\\gamma^0_X/\gamma^0_V-1}', '\\gamma^1_X','{\\gamma^1_X/\\gamma^1_V-1}',\
   '\\delta_X','{{\\delta_X/\\delta_I-1}}','\\sigma_X']
sigfig = [4,1,3,2,4,2,2,2,2,2,2,2,3]
for p,pn, s in zip(pars,pars_n,sigfig):
    print '${}$'.format(pn)
    for i in xrange(5):
        print '&'
        if pn[0] == 'R':
            dum = numpy.percentile(fit[p][:,i]/(fit[p][:,1]-fit[p][:,2])[:,None],(50,50-34,50+34)) 
        elif pn[0:2] == '{{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,4]-1,(50,50-34,50+34)) 
        elif pn[0] == '{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34)) 
        else :
            dum = numpy.percentile(fit[p][:,i],(50,50-34,50+34))            
        print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'

filts = ['U','B','V','R','I']
for i in xrange(5):
    for j in xrange(i+1,5):
        print '{}-{} {} {}'.format(filts[i],filts[j], \
            numpy.std((fit['gamma'][:,i]-fit['gamma'][:,j])[:,None]*fit['k']),numpy.std((fit['rho1'][:,i]-fit['rho1'][:,j])[:,None]*fit['R']))


mega = numpy.array([fit['Delta'],fit['EW'][:,:,0],fit['EW'][:,:,1],fit['sivel'], \
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),((fit['gamma1'][:,1] - fit['gamma1'][:,2])[:,None]*fit['k1']),(fit['rho1'][:,4][:,None]*fit['R'])])

# dum = numpy.corrcoef(mega)
# print "observable correlation coefficients"
# print " \\\\\n".join([" & ".join(map('{0:.2f}'.format, line)) for line in dum])



# mega = numpy.array([fit['Delta'],fit['EW'][:,:,0],fit['EW'][:,:,1],fit['sivel'], \
#     ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R'])])

# mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), \
#     ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']).flatten(),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']).flatten()])

corrarray = []

for i in xrange(mega.shape[1]):
    corrarray.append(numpy.corrcoef(mega[:,i,:]))    

corrarray=numpy.array(corrarray)

dum1, dumm, dump = numpy.percentile(corrarray,(50,50-34,50+34),axis=0)

# dum = numpy.zeros()
# dum = numpy.corrcoef(mega)
print "observable correlation coefficients"
# dum=numpy.zeros((6,18))
for i1 in xrange(7):
    for i2 in xrange(7):
        print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1[i1,i2],dump[i1,i2]-dum1[i1,i2],dumm[i1,i2]-dum1[i1,i2]),
        if (i2 != 6):
            print "&",

    print "\\\\" 

wefwe




# pars = ['gamma','rho1']
# pars_n = ['R_{','R_{\\delta ']
# sigfig = [3,3]
# for p,pn, s in zip(pars,pars_n,sigfig):
#     print '${}X}}$'.format(pn)
#     for i in xrange(5):
#         dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34))            
#         print  '${1:6.{0}f}^{{{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
#     print '\\\\'


ebvdelta  = (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']
(y,ymin,ymax) = numpy.percentile(ebvdelta,(50,50-34,50+34),axis=0)
wmax = numpy.argmax(y)
wmin = numpy.argmin(y)
ydelta = y
print 'Extreme values of E_delta(B-V)'

print "{:6.2f}_{{{:6.2f}}}^{{+{:6.2f}}}".format(y[wmax],ymin[wmax]-y[wmax],ymax[wmax]-y[wmax])
print "{:6.2f}_{{{:6.2f}}}^{{+{:6.2f}}}".format(y[wmin],ymin[wmin]-y[wmin],ymax[wmin]-y[wmin])

ebvdelta  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']
(y,ymin,ymax) = numpy.percentile(ebvdelta,(50,50-34,50+34),axis=0)
wmax = numpy.argmax(y)
wmin = numpy.argmin(y)
print 'Extreme values of E_gamma(B-V)'
print "{:6.2f}_{{{:6.2f}}}^{{+{:6.2f}}}".format(y[wmax],ymin[wmax]-y[wmax],ymax[wmax]-y[wmax])
print "{:6.2f}_{{{:6.2f}}}^{{+{:6.2f}}}".format(y[wmin],ymin[wmin]-y[wmin],ymax[wmin]-y[wmin])



pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()


sivel,sivel_err,x1,x1_err,zcmb = sivel.sivel(data)

use = numpy.isfinite(sivel)

# #  The ordering is 'Ca','Si','U','B','V','R','I'

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


# ebvdelta  = (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']
# ebvgamma  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']

# for s,a,b in zip(snname, numpy.median(ebvdelta,axis=0),numpy.median(ebvgamma,axis=0)):
#   print s,a,b



# for i in xrange(len(sivel)):
#     print '{0} & ${1:5.1f} \pm {2:3.1f}$ & ${3:5.1f} \pm {4:3.1f}$& ${7:5.0f} \pm {8:3.0f}$ & ${5[0]:6.2f} \pm {6[0]:6.2f}$ & ${5[1]:6.2f} \pm {6[1]:6.2f}$& ${5[2]:6.2f} \pm {6[2]:6.2f}$& ${5[3]:6.2f} \pm {6[3]:6.2f}$& ${5[4]:6.2f} \pm {6[4]:6.2f}$ \\\\'.format(snname[i], EW_obs[i,0], numpy.sqrt(EW_cov[i,0,0]),
#         EW_obs[i,1], numpy.sqrt(EW_cov[i,1,1]), mag_obs[i,:], numpy.sqrt(numpy.diagonal(mag_cov[i,:,:])),sivel[i],sivel_err[i])

