#!/usr/bin/env python

import pickle
import cPickle
import numpy
import sivel
import sncosmo
import fitz_band

f = open('fix3_x1_decorr.pkl','rb')
# (fit, _) = pickle.load(f)
fit = pickle.load(f)
f.close()


for key in fit.keys():
    print key, fit[key].min(), fit[key].max()

#special variables

# ev flipped
fixev = fit['ev']
fixev = -fixev * numpy.sign(fixev[:,2])[:,None]


print "projection of ev out of gamma plane"
proj = []
for ev, evorig, g1, g2 in zip(fixev, fit['ev'], fit['gamma'], fit['rho1']):
    y=numpy.array([numpy.dot(evorig,g1), numpy.dot(evorig,g2)])
    M = numpy.array([[numpy.dot(g1,g1), numpy.dot(g1,g2)], \
        [numpy.dot(g2,g1), numpy.dot(g2,g2)]])
    com = numpy.dot(numpy.linalg.inv(M),y)
    proj.append(1-numpy.linalg.norm(com[0]*g1+com[1]*g2)**2 / numpy.linalg.norm(evorig)**2)
dum1, dumm, dump =  numpy.percentile(proj,(50,50-34,50+34),axis=0)
print "${:.3f}^{{+{:.3f}}}_{{{:.3f}}}$".format(dum1,dump-dum1,dumm-dum1)

print "ev_sig size"
dum1, dumm, dump =  numpy.percentile(fit['ev_sig'],(50,50-34,50+34),axis=0)
print "{:.3f}^{{+{:.3f}}}_{{{:.3f}}}".format(dum1,dump-dum1,dumm-dum1)

print "global parameter covariance matrix"
nlinks = fit['gamma'].shape[0]
mega = [] 
for i in xrange(5):
    mega.append([fit['alpha'][:,i],fit['beta'][:,i],fit['eta'][:,i],\
        fit['zeta'][:,i],fit['gamma'][:,i],fit['rho1'][:,i],fit['ev_sig'][:]*fit['ev'][:,i]])
mega = numpy.array(mega)
mega = numpy.reshape(mega, (35,nlinks))
print mega.shape
ans = numpy.cov(mega)
def frexp10(x):
    exp = int(numpy.floor(numpy.log10(abs(x))))
    return x / 10**exp, exp

for i1 in xrange(35):
    for i2 in xrange(35):
        print "{0}^{{{1:+03}}}".format(*frexp10(ans[i1,i2])),
        if (i2 != 34):
            print "&",
    print "\\\\" 
wefew



print "the table"

pars = ['alpha','alpha','beta','beta','eta','eta','zeta','zeta','gamma','gamma','rho1','rho1']
pars_n = ['\\alpha_X','{\\alpha_X}/\\alpha_{\hat{V}}-1','\\beta_X','{\\beta_X}/\\beta_{\hat{V}}-1',\
  '\\eta_X','{\\eta_X}/\\eta_{\\hat{V}}-1', '\\zeta_X','{\\zeta_X}/\\zeta_{\\hat{V}}-1','\\gamma^0_X', '{\\gamma^0_X}/\gamma^0_{\\hat{V}}-1', '\\gamma^1_X','{\\gamma^1_X}/\\gamma^1_{\\hat{V}}-1']
sigfig = [4,1,3,2,4,2,2,2,2,2,2,2,3]
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

# the smallest ev_sig
print "ev_sig"
dum = numpy.percentile(fit['ev_sig'],(50,50-34,50+34))
print  '${1:6.{0}f}^{{+{2:6.{0}f}}}_{{{3:6.{0}f}}}$'.format(3,dum[0], dum[2]-dum[0],dum[1]-dum[0] )

print "ev_sig min"
print  fit['ev_sig'].min()

print "standard deviation of delta"
print fit['Delta'].std()



print "standard deviations of E(B-V)"
print numpy.std(fit['k']*((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None])
print numpy.std(fit['R']*((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None])

print numpy.std((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'])
print numpy.std((fit['ev_sig']*(fit['ev'][:,1]-fit['ev'][:,2]))[:,None]*fit['mag_int_raw'])

filts = ['U','B','V','R','I']
for i in xrange(5):
    for j in xrange(i+1,5):
        print '{}-{} {} {}'.format(filts[i],filts[j], \
            numpy.std((fit['gamma'][:,i]-fit['gamma'][:,j])[:,None]*fit['k']),numpy.std((fit['rho1'][:,i]-fit['rho1'][:,j])[:,None]*fit['R']))


mega = numpy.array([fit['Delta'],fit['EW'][:,:,0],fit['EW'][:,:,1],fit['sivel'], fit['x1'],\
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']), fit['mag_int_raw']*(fit['ev_sig']*fit['ev'][:,2])[:,None]])

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
for i1 in xrange(8):
    for i2 in xrange(8):
        print "{:.2f}^{{+{:.2f}}}_{{{:.2f}}}".format(dum1[i1,i2],dump[i1,i2]-dum1[i1,i2],dumm[i1,i2]-dum1[i1,i2]),
        if (i2 != 7):
            print "&",

    print "\\\\" 


# Effecrive RB and other numbers
print "Effecrive RB and other numbers"
print 'RB = AB / (AB-AV)'
dum1, dumm, dump = numpy.percentile(1/(1-fit['ev'][:,2]/fit['ev'][:,1]),(50,50-34,50+34))
print "{:.1f}^{{+{:.1f}}}_{{{:.1f}}}".format(dum1,dump-dum1,dumm-dum1)

print 'RB = AB / (AB-AR)'
dum1, dumm, dump = numpy.percentile(1/(1-fit['ev'][:,3]/fit['ev'][:,1]),(50,50-34,50+34))
print "{:.1f}^{{+{:.1f}}}_{{{:.1f}}}".format(dum1,dump-dum1,dumm-dum1)

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

#ratio phi_b/(phi_b - phi_v)
print 'ratio phi_b/(phi_b - phi_v)'
(y,ymin,ymax) = numpy.percentile(fit['ev'][:,1]/(fit['ev'][:,1]-fit['ev'][:,2]),(50,50-34,50+34),axis=0)
print "{:6.1f}_{{{:6.1f}}}^{{+{:6.1f}}}".format(y,ymin-y,ymax-y)

print 'ratio phi_b/(phi_b - phi_r)'
(y,ymin,ymax) = numpy.percentile(fit['ev'][:,1]/(fit['ev'][:,1]-fit['ev'][:,3]),(50,50-34,50+34),axis=0)
print "{:6.1f}_{{{:6.1f}}}^{{+{:6.1f}}}".format(y,ymin-y,ymax-y)


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel, sivel_err,x1,x1_err,_,_,_,EWFe4800 = sivel.sivel(data)


use = numpy.isfinite(sivel)

# #  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
x1 = x1[use]
x1_err = x1_err[use]
EWFe4800=EWFe4800[use]
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]
nsne= len(use)
shit = []
for m in mag_cov:
    shit.append(numpy.sqrt(numpy.diag(m)))
print numpy.median(shit)

shit=[]
for i in xrange(fit['mag_int_raw'].shape[0]):
    shit.append(numpy.corrcoef(numpy.array([EWFe4800, fit['mag_int_raw'][i,:]]))[0,1])
(y,ymin,ymax) = numpy.percentile(shit,(50,50-34,50+34),axis=0)
print y, y-ymin,ymax-y

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
shit = []
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))
    shit.append(numpy.sqrt(numpy.diag(color_cov[i])))

print numpy.max(shit)
print numpy.median(shit)

# ebvdelta  = (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None] * fit['R']
# ebvgamma  = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None] * fit['k']

# for s,a,b in zip(snname, numpy.median(ebvdelta,axis=0),numpy.median(ebvgamma,axis=0)):
#   print s,a,b


snname=numpy.array(data['snlist'])[use]

print 'outputs'

EW_mn = EW_obs.mean(axis=0)
sivel_mn = sivel.mean()
for i in xrange(len(sivel)):
    (y0,ymin0,ymax0) = numpy.percentile(fit['EW'][:,i,0],(50,50-34,50+34),axis=0)
    (y1,ymin1,ymax1) = numpy.percentile(fit['EW'][:,i,1],(50,50-34,50+34),axis=0)
    (y2,ymin2,ymax2) = numpy.percentile(fit['sivel'][:,i],(50,50-34,50+34),axis=0)
    (y3,ymin3,ymax3) = numpy.percentile(fit['x1'][:,i],(50,50-34,50+34),axis=0)
    (y4,ymin4,ymax4) = numpy.percentile(fit['Delta'][:,i],(50,50-34,50+34),axis=0)
    (y5,ymin5,ymax5) = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])*fit['k'][:,i],(50,50-34,50+34),axis=0)
    (y6,ymin6,ymax6) = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])*fit['R'][:,i],(50,50-34,50+34),axis=0)
    (y7,ymin7,ymax7) = numpy.percentile(fit['ev_sig']*fit['ev'][:,2]*fit['mag_int_raw'][:,i],(50,50-34,50+34),axis=0)
    print '{0} & ${1:5.1f}^{{+{2:3.1f}}}_{{-{3:3.1f}}}$ & ${4:5.1f}^{{+{5:3.1f}}}_{{-{6:3.1f}}}$ & ${7:5.0f}^{{+{8:3.0f}}}_{{-{9:3.0f}}}$ & ${10:6.2f}^{{+{11:6.2f}}}_{{-{12:6.2f}}}$ & ${13:5.3f}^{{+{14:5.3f}}}_{{-{15:5.3f}}}$  & ${16:5.3f}^{{+{17:5.3f}}}_{{-{18:5.3f}}}$ & ${19:5.3f}^{{+{20:5.3f}}}_{{-{21:5.3f}}}$ & ${19:5.3f}^{{+{20:5.3f}}}_{{-{21:5.3f}}}$\\\\'.format(\
        snname[i], EW_mn[0]+y0,ymax0-y0, y0-ymin0, EW_mn[1]+y1,ymax1-y1, y1-ymin1,
        sivel_mn+ y2, ymax2-y2, y2-ymin2,y3, ymax3-y3, y3-ymin3,y4, ymax4-y4, y4-ymin4,
        y5, ymax5-y5, y5-ymin5, y6, ymax6-y6, y6-ymin6, y7, ymax7-y7, y7-ymin7
        )
    #, numpy.sqrt(EW_cov[i,0,0]),
     #   EW_obs[i,1], numpy.sqrt(EW_cov[i,1,1]), mag_obs[i,:], numpy.sqrt(numpy.diagonal(mag_cov[i,:,:])),sivel[i],sivel_err[i], x1[i], x1_err[i])
# ${5[0]:6.2f} \pm {6[0]:6.2f}$ & ${5[1]:6.2f} \pm {6[1]:6.2f}$& ${5[2]:6.2f} \pm {6[2]:6.2f}$& ${5[3]:6.2f} \pm {6[3]:6.2f}$& ${5[4]:6.2f} \pm {6[4]:6.2f}$ & ${9:6.2f} \pm {10:6.2f}$
werfwe

print 'inputs'
for i in xrange(len(sivel)):
    print '{0} & ${1:5.1f} \pm {2:3.1f}$ & ${3:5.1f} \pm {4:3.1f}$& ${7:5.0f} \pm {8:3.0f}$ & ${5[0]:6.2f} \pm {6[0]:6.2f}$ & ${5[1]:6.2f} \pm {6[1]:6.2f}$& ${5[2]:6.2f} \pm {6[2]:6.2f}$& ${5[3]:6.2f} \pm {6[3]:6.2f}$& ${5[4]:6.2f} \pm {6[4]:6.2f}$ & ${9:6.2f} \pm {10:6.2f}$\\\\'.format(snname[i], EW_obs[i,0], numpy.sqrt(EW_cov[i,0,0]),
        EW_obs[i,1], numpy.sqrt(EW_cov[i,1,1]), mag_obs[i,:], numpy.sqrt(numpy.diagonal(mag_cov[i,:,:])),sivel[i],sivel_err[i], x1[i], x1_err[i])

