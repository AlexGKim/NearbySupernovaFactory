#!/usr/bin/env python

import pickle
import cPickle
import numpy
import sivel
import sncosmo
import fitz_band

f = open('fix3_decorr.pkl','rb')
fit = pickle.load(f)
f.close()

for key in fit.keys():
    print key, fit[key].min(), fit[key].max()


#special variables

# ev flipped
fixev = fit['ev']
fixev = fixev * numpy.sign(fixev[:,4])[:,None]


print "projection of ev outside of gamma plane"
proj = []
for ev, evorig, g1, g2 in zip(fixev, fit['ev'], fit['gamma'], fit['rho1']):
    y=numpy.array([numpy.dot(evorig,g1), numpy.dot(evorig,g2)])
    M = numpy.array([[numpy.dot(g1,g1), numpy.dot(g1,g2)], \
        [numpy.dot(g2,g1), numpy.dot(g2,g2)]])
    com = numpy.dot(numpy.linalg.inv(M),y)
    proj.append(1-numpy.linalg.norm(com[0]*g1+com[1]*g2)**2 / numpy.linalg.norm(evorig)**2)
dum1, dumm, dump =  numpy.percentile(proj,(50,50-34,50+34),axis=0)
print "${:.2f}^{{+{:.2f}}}_{{{:.2f}}}$".format(dum1,dump-dum1,dumm-dum1)

print "ev_sig size"
dum1, dumm, dump =  numpy.percentile(fit['ev_sig'],(50,50-34,50+34),axis=0)
print "{:.3f}^{{+{:.3f}}}_{{{:.3f}}}".format(dum1,dump-dum1,dumm-dum1)

print "the table"

pars = ['alpha','alpha','beta','beta','eta','eta','gamma','gamma','rho1','rho1']
pars_n = ['\\alpha_X','{\\alpha_X}/\\alpha_{\\hat{V}-1}','\\beta_X','{\\beta_X}/\\beta_{\\hat{V}-1}',\
  '\\eta_X','{\\eta_X}/\\eta_{\\hat{V}-1}', '\\gamma^0_X', '{\\gamma^0_X}/\gamma^0_{\\hat{V}-1}', '\\gamma^1_X','{\\gamma^1_X}/\\gamma^1_{\\hat{V}-1}']
sigfig = [4,1,3,2,4,2,2,2,2,2,3]
for p,pn, s in zip(pars,pars_n,sigfig):
    print '${}$'.format(pn)
    for i in xrange(5):
        print '&',
        if pn[0] == 'R':
            dum = numpy.percentile(fit[p][:,i]/(fit[p][:,1]-fit[p][:,2])[:,None],(50,50-34,50+34)) 
        elif pn[0] == '{':
            dum = numpy.percentile(fit[p][:,i]/fit[p][:,2]-1,(50,50-34,50+34)) 
        else :
            dum = numpy.percentile(fit[p][:,i],(50,50-34,50+34))            
        print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(s,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
    print '\\\\'

print '$\sigma_p \phi_X$'
for i in xrange(5):
    print '&',
    dum = numpy.percentile(fit['ev_sig']*fixev[:,i],(50,50-34,50+34))
    print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(3,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
print '\\\\'
print '${\\phi_X/\\phi_{\\hat{V}}-1}$'
for i in xrange(5):
    print '&',
    dum = numpy.percentile(fit['ev'][:,i]/fit['ev'][:,2]-1,(50,50-34,50+34))
    print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(3,dum[0], dum[2]-dum[0],dum[1]-dum[0] )
print '\\\\'



print "standard deviation of delta"
print (fit['Delta']-fit['Delta'][:,0][:,None])[:,1:].shape
print (fit['Delta']-fit['Delta'][:,0][:,None])[:,1:].std()
print fit['Delta'].std()

# trans = [[0,0,1.,0,0],[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
# trans = numpy.array(trans)
# si=[]
# cmat = []
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     covmat = numpy.dot(trans,numpy.dot(numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T),trans.T))
#     D = numpy.sqrt(numpy.diag(covmat))

#     R = covmat/numpy.outer(D,D)

#     si.append(D)
#     cmat.append(R)

# si = numpy.array(si)
# cmat = numpy.array(cmat)


# dum1, dumm, dump =  numpy.percentile(si,(50,50-34,50+34),axis=0)
# for i2 in xrange(5):
#     print "{:.3f}^{{+{:.3f}}}_{{{:.3f}}}".format(dum1[i2],dump[i2]-dum1[i2],dumm[i2]-dum1[i2])



# dum1, dumm, dump = numpy.percentile(cmat,(50,50-34,50+34),axis=0)

# # dum = numpy.zeros()
# # dum = numpy.corrcoef(mega)
# print "intrinsic correlation coefficients"
# # dum=numpy.zeros((6,18))
# for i1 in xrange(5):
#     for i2 in xrange(5):
#         print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1[i1,i2],dump[i1,i2]-dum1[i1,i2],dumm[i1,i2]-dum1[i1,i2]),
#         if (i2 != 4):
#             print "&",

#     print "\\\\" 








# dum = numpy.zeros((5,5))
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

# dum/= len(fit['L_Omega'])
# print "average L_Omega"
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

# dum = numpy.zeros((5,5))
# for x1 in fit['L_Omega']:
#     dum= dum+ numpy.dot(x1,x1.T)

# dum/= len(fit['L_Omega'])
# print "average correlation"
# print " \\\\\n".join([" & ".join(map('{0:.2f}'.format, line)) for line in dum])

# trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
# trans = numpy.array(trans)

# dum = numpy.zeros((4,4))
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     dum= dum+ numpy.dot(trans,numpy.dot(numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T),trans.T))
# dum/= len(fit['L_Omega'])

# # color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
# # print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

# print "color covariance"
# dumsig = numpy.sqrt(numpy.diag(dum))
# print [" , ".join(map('{0:.3f}'.format, dumsig))]
# dumcor =  dum/ numpy.dot(dumsig[:,None],dumsig[None,:])
# print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in dumcor])

print "standard deviations of E(B-V)"
print numpy.std(fit['k']*((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None])
print numpy.std(fit['R']*((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None])



filts = ['U','B','V','R','I']
for i in xrange(5):
    for j in xrange(i+1,5):
        print '{}-{} {} {}'.format(filts[i],filts[j], \
            numpy.std((fit['gamma'][:,i]-fit['gamma'][:,j])[:,None]*fit['k']),numpy.std((fit['rho1'][:,i]-fit['rho1'][:,j])[:,None]*fit['R']))


mega = numpy.array([fit['Delta'],fit['EW'][:,:,0],fit['EW'][:,:,1],fit['sivel'], \
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R'])])

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
for i1 in xrange(6):
    for i2 in xrange(6):
        print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1[i1,i2],dump[i1,i2]-dum1[i1,i2],dumm[i1,i2]-dum1[i1,i2]),
        if (i2 != 5):
            print "&",

    print "\\\\" 

#         # dum[i1,3*i2]=dum1[i1,i2]
#         # dum[i1,3*i2+1]=dump[i1,i2]-dum1[i1,i2]
#         # dum[i1,3*i2+2]=dum1[i1,i2]-dumm[i1,i2]



# # dum = zip(dum1,dump-dum1,dum1-dumm)
# # dum = numpy.reshape(numpy.array(dum),(6,18))
# print " \\\\\n".join([" & ".join(map('{:.2f}^{:.2f}_{:.2f}'.format,line)) for line in dum])


# print " \\\\\n".join([" & ".join(map('${0:.2f}^1_1$'.format,line)) for line in dum])



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

sivel, sivel_err,x1,x1_err,_,_,_,EWFe4800,_ = sivel.sivel(data)

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


snname=numpy.array(data['snlist'])[use]

for i in xrange(len(sivel)):
    print '{0} & ${1:5.1f} \pm {2:3.1f}$ & ${3:5.1f} \pm {4:3.1f}$& ${7:5.0f} \pm {8:3.0f}$ & ${5[0]:6.2f} \pm {6[0]:6.2f}$ & ${5[1]:6.2f} \pm {6[1]:6.2f}$& ${5[2]:6.2f} \pm {6[2]:6.2f}$& ${5[3]:6.2f} \pm {6[3]:6.2f}$& ${5[4]:6.2f} \pm {6[4]:6.2f}$ \\\\'.format(snname[i], EW_obs[i,0], numpy.sqrt(EW_cov[i,0,0]),
        EW_obs[i,1], numpy.sqrt(EW_cov[i,1,1]), mag_obs[i,:], numpy.sqrt(numpy.diagonal(mag_cov[i,:,:])),sivel[i],sivel_err[i])



import json, codecs
EW_mn = EW_obs.mean(axis=0)
sivel_mn = sivel.mean()
def convert(fit,EW_mn,sivel_mn):
    rm = ['c_raw','alpha_raw','beta_raw','eta_raw','gamma01','gamma02','gamma03','gamma04','gamma05','k_unit','lp__', \
        'Delta_scale','Delta_unit','R_unit','rho11','rho12','rho13','rho14','rho15','EW','sivel','gamma','rho1','k','R', \
        'ev_sig','ev','mag_int_raw','mag_int']
    use = dict(fit)

    use['EWCa'] = use['EW'][:,:,0]+EW_mn[0]
    use['EWSi'] = use['EW'][:,:,1]+EW_mn[1]
    use['lSi']=use['sivel'] + sivel_mn
    use['gamma0']=use['gamma']
    use['gamma1']=use['rho1']
    use['g0'] = use['k']
    use['g1'] = use['R']
    use['sigmap']=use['ev_sig']
    use['phi']=use['ev']
    use['p'] = use['mag_int_raw']
    for r in rm:
        del use[r]

    for key in use.keys():
        use[key] = use[key].tolist()

    json.dump(use, codecs.open('mIIchain.json', 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4) ### this saves the array in .json format
convert(fit, EW_mn, sivel_mn)
wefwe