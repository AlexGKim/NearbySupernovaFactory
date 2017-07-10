#!/usr/bin/env python

import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo
import scipy
import cPickle
import matplotlib as mpl
import sivel
import flip

dirname = 'output25/'

mpl.rcParams['font.size'] = 18

f = open('temp25.pkl','rb')
(fit,_) = pickle.load(f)

# flipin=[]
# for i in xrange(8):
#     flipin.append(numpy.median(fit['rho1'][i*2500:(i+1)*2500,0]))

# flipin = numpy.array(flipin)

# for key in fit.keys():
#     shit = numpy.percentile(fit[key],(50,50-34,50+34),axis=0) 
#     print key, shit

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel,sivel_err,x1,x1_err,zcmb,zerr = sivel.sivel(data)


use = numpy.isfinite(sivel)

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

sivel=sivel[use]
sivel_err = sivel_err[use]
x1=x1[use]
x1_err = x1_err[use]
zcmb = zcmb[use]
zerr = zerr[use]

EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

snnames = numpy.array(data['snlist'])[use]

EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn

nsne, nmags = mag_obs.shape
color_obs = numpy.zeros((nsne,nmags-1))
color_obs[:,0] = mag_renorm[:,0]- mag_renorm[:,2]
color_obs[:,1] = mag_renorm[:,1]- mag_renorm[:,2]
color_obs[:,2] = mag_renorm[:,3]- mag_renorm[:,2]
color_obs[:,3] = mag_renorm[:,4]- mag_renorm[:,2]

EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

pkl_file.close()

plt.plot(fit['rho1'][::10,:])
plt.title(r'$\delta')
pp = PdfPages(dirname+"/delta_chain.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))

(y, ymin, ymax) = numpy.percentile(fit['Delta'],(50,50-34,50+34),axis=0)

import rdata
rd = rdata.rdata()
x=[]
dx=[]
for dum in snnames:
  if dum in rd:
    x.append(rd[dum][0])
    dx.append(rd[dum][1])
  else:
    x.append(numpy.nan)
    dx.append(numpy.nan)
x=numpy.array(x)
dx=numpy.array(dx)
w = numpy.isfinite(x)
plt.errorbar(x[w], y[w], xerr=[dx[w],dx[w]],yerr=[y[w]-ymin[w],ymax[w]-y[w]],fmt='o')
plt.xlabel(r'Hubble Residual')
plt.ylabel(r'$\Delta$')
pp = PdfPages(dirname+'/DeltamB.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.errorbar(zcmb, y, yerr=[y-ymin,ymax-y],fmt='o')
plt.xlabel(r'$z_{cmb}$')
plt.ylabel(r'$\Delta$')
pp = PdfPages(dirname+'/Deltaz.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

dmu = numpy.sqrt((5./numpy.log(10)*300./(zcmb*3e5))**2 + ((ymax-ymin)/2)**2)
smallz = zerr < 0.001

pull = y/dmu

plt.hist([pull[smallz],pull[numpy.logical_not(smallz)]],log=True,normed=False,bins=20,label=[r'low $\sigma_z$',r'not low $\sigma_z$'])
# print snnames[pull > 2]
plt.legend()

plt.xlabel(r'$\Delta / \sqrt{(\sigma^2_{\Delta} + \sigma^2_{\Delta v})}$')
pp = PdfPages(dirname+'/Deltanorm_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(x, xmin, xmax) = numpy.percentile( ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),(50,50-34,50+34),axis=0)
(y, ymin, ymax) = numpy.percentile(((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']),(50,50-34,50+34),axis=0)

plt.errorbar(x,y,xerr=(x-xmin,xmax-x), yerr=(y-ymin,ymax-y),fmt='o')
plt.xlim((-0.08,0.05))
# plt.ylim((-0.02,0.05))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$E_\delta(B-V)$')
pp = PdfPages(dirname+"/egammaedelta_corner.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(y, ymin, ymax) = numpy.percentile(fit['EW'][:,:,1],(50,50-34,50+34),axis=0)
plt.errorbar(x1, y, xerr=[x1_err,x1_err],yerr=[y-ymin,ymax-ymin],fmt='o')
plt.xlabel(r'$x_1$')
plt.ylabel(r'$EW_{Si}$')
plt.tight_layout()
pp = PdfPages(dirname+"/x1si.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(y, ymin, ymax) = numpy.percentile(fit['rho1'][:,0][:,None]*fit['R'],(50,50-34,50+34),axis=0)
plt.errorbar(x1, y, xerr=[x1_err,x1_err],yerr=[y-ymin,ymax-ymin],fmt='o',alpha=0.2)
plt.scatter(x1, y,s=1,alpha=0.8)
plt.xlabel(r'$x_1$')
plt.ylabel(r'$A_{\delta U}$')
plt.tight_layout()
pp = PdfPages(dirname+"/x1D.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# offset = 0.01
# au = fit['gamma'][:,0][:,None]*fit['k'] + fit['rho1'][:,0][:,None]*fit['R']
# ab = fit['gamma'][:,1][:,None]*fit['k'] + fit['rho1'][:,1][:,None]*fit['R']
# av = fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R']


# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0)
# (y, ymin, ymax) = numpy.percentile(au-av,(50,50-34,50+34),axis=0)

# plt.scatter(x ,numpy.median(av/(ab-av),axis=0))
# plt.ylim((-2,6))
# plt.xlim((-0.1,.4))
# plt.ylabel(r'$R_{eff,V}$')
# plt.xlabel(r'$E_{eff}(B-V)$')
# pp = PdfPages(dirname+"/Rveff4.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# x=x-x.min()-offset*1.3
# y=y-y.min()-offset
# plt.scatter(x,y)
# # x_=numpy.array([-0.03,0.45])
# # plt.plot(x_,(4.87-2.96)*x_)
# plt.xlim((-0.04,0.5))
# plt.ylim((-0.04,1))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$E_{eff}(U-V)$')
# pp = PdfPages(dirname+"/Rveff.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (y, ymin, ymax) = numpy.percentile((au-av-y.min()-offset)/(ab-av-x.min()-offset),(50,50-34,50+34),axis=0)
# plt.scatter(x,y)
# # plt.ylim((1.8,2.3))
# #plt.xlim((0,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$R_{eff,U} - R_{eff,V}$')
# pp = PdfPages(dirname+"/Rveff2.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# temp = (fit['gamma'][:,1][:,None]*fit['k']-fit['gamma'][:,2][:,None]*fit['k'] ) /\
#  (fit['rho1'][:,1][:,None]*fit['R'] - fit['rho1'][:,2][:,None]*fit['R'])
# (x, xmin, xmax) = numpy.percentile(temp,(50,50-34,50+34),axis=0)
# plt.scatter(x,y)
# # plt.ylim((1.8,2.3))
# plt.xlabel(r'$E_\gamma(B-V)/E_\rho(B-V)$')
# plt.ylabel(r'$R_{eff,U} - R_{eff,V}$')
# pp = PdfPages(dirname+"/Rveff3.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


correction = [fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
    + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']\
    + fit['gamma'][:,i][:,None]*fit['k'] + fit['gamma1'][:,i][:,None]*fit['k1'] + fit['rho1'][:,i][:,None]*fit['R']\
    for i in xrange(5)]
correction = numpy.array(correction)
correction = correction - correction[2,:,:]
# correction_median = numpy.median(correction,axis=1)
from matplotlib.ticker import NullFormatter
cind=[0,1,3,4]
cname = ['U','B','R','I']
mpl.rcParams['font.size'] = 14

fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'width_ratios':[1,.2]})
for i in xrange(4):
    (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
    err = numpy.sqrt(color_cov[:,i,i] + ((ymax-ymin)/2)**2)
    axes[i,0].errorbar(color_obs[:,i],color_obs[:,i]-y,xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[err,err],fmt='.',alpha=0.4)

    # axes[i,0].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]-y,xerr=[y-ymin,ymax-y], yerr=[err,err],fmt='.',alpha=0.4)
    # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.5)

 
    # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i],xerr=[y-ymin,ymax-y], yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])],fmt='.',alpha=0.4)
    miny = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).min()
    maxy = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).max()
    # axes[i].plot([miny,maxy],[miny,maxy])
    axes[i,0].set_xlabel(r'$({0}_o-V_o)$'.format(cname[i]))
    lname = r'$\Delta({0}-V)$'.format(cname[i])
    axes[i,0].set_ylabel(lname)
    axes[i,1].hist(color_obs[:,i]-y, orientation='horizontal')
    axes[i,1].set_ylim(axes[i,0].get_ylim())
    axes[i,1].yaxis.set_major_formatter(NullFormatter())
    axes[i,1].xaxis.set_major_formatter(NullFormatter())
    # axes[i].axhline(y=0,linestyle=':')
fig.subplots_adjust(hspace=.4, wspace=.04)
fig.set_size_inches(8,11)
filename = dirname+'/residual.pdf'
pp = PdfPages(filename)
plt.savefig(pp,format='pdf')
pp.close()



# correction_gerard = [fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
#     + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']\
#     + fit['gamma'][:,i][:,None]*fit['k']+ fit['gamma1'][:,i][:,None]*fit['k1'] + fit['rho1'][:,i][:,None]*fit['R']\
#     for i in xrange(5)]

# correction_gerard = numpy.array(correction_gerard)
# correction_gerard = correction_gerard - correction_gerard[1,:,:]

# gin =3
# gerard_obs = numpy.zeros((nsne,nmags-1))
# gerard_obs[:,0] = mag_renorm[:,0]- mag_renorm[:,1]
# gerard_obs[:,1] = mag_renorm[:,2]- mag_renorm[:,1]
# gerard_obs[:,2] = mag_renorm[:,3]- mag_renorm[:,1]
# gerard_obs[:,3] = mag_renorm[:,4]- mag_renorm[:,1]

# trans = [[1.,-1.,0,0,0],[0.,1,0,-1,0],[0.,0,1,-1,0],[0.,0,0,-1,0]]
# trans = numpy.array(trans)
# gerard_cov = numpy.zeros((nsne,4,4))
# for i in xrange(nsne):
#     gerard_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))

# cind=[0,1,3,4]
# cname = ['U','B','R','I']
# fig, axes = plt.subplots(nrows=4)
# (x, xmin, xmax) = numpy.percentile(correction_gerard[0,:,:],(50,50-34,50+34),axis=0)
# xerr = numpy.sqrt(gerard_cov[:,0,0] + ((ymax-ymin)/2)**2)   
# for i in xrange(4):
#     (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
#     yerr = numpy.sqrt(gerard_cov[:,i,i] + ((ymax-ymin)/2)**2)
#     axes[i].errorbar(gerard_obs[:,0]-x,color_obs[:,i]-y,xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.',alpha=0.5)
#     axes[i].set_xlabel(r'$({0}_o-B_o) - ({0}-B)_{{model}}$'.format('U'))
#     lname = r'$({0}_o-V_o) - ({0}-V)_{{model}}$'.format(cname[i])
#     axes[i].set_ylim((-0.2,0.2))
#     axes[i].set_ylabel(lname)
#     axes[i].axhline(y=0,linestyle=':')
# fig.subplots_adjust(hspace=.3)
# fig.set_size_inches(8,11)
# filename = dirname+'/residual_gerard.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()

# wefwe


cind=[0,1,3,4]
cname = ['U','B','R','I']
fig, axes = plt.subplots(nrows=4,sharex=True)
for i in xrange(4):
    (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
    err = numpy.sqrt(color_cov[:,i,i] + ((ymax-ymin)/2)**2)
    axes[i].errorbar(x1,color_obs[:,i]-y,xerr=[x1_err,x1_err], yerr=[err,err],fmt='.',alpha=0.2)
    axes[i].scatter(x1,color_obs[:,i]-y,s=1,alpha=0.8)
    lname = r'$\Delta ({0}-V)$'.format(cname[i])
    axes[i].set_ylabel(lname)
    axes[i].axhline(y=0,linestyle=':')

axes[3].set_xlabel(r'$x_1$')

# fig.subplots_adjust(hspace=.3)
fig.set_size_inches(10,11)
plt.setp(axes[0].get_xticklabels(),visible=False)
plt.setp(axes[1].get_xticklabels(),visible=False)
plt.setp(axes[2].get_xticklabels(),visible=False)
plt.subplots_adjust(hspace=.1)
#plt.tight_layout()
filename = dirname+'/residualx1.pdf'
pp = PdfPages(filename)
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
plt.clf()

# plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0), mag_obs[:,2])
# x=numpy.array([-0.08, 0.37])
# plt.plot(x,-28.6+x*3.96)
# plt.show()


# dum = numpy.zeros((5,5))
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

# dum/= len(fit['L_Omega'])
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

# trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
# trans = numpy.array(trans)
# color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in color_cov])

# correction = [ fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
#     + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] \
#     + fit['gamma'][:,i][:,None]*fit['k'] \
#     + fit['Delta'] \
#     for i in xrange(5)]

# correction = numpy.array(correction)
# correction_median = numpy.median(correction,axis=1)

# outlier  = numpy.where((mag_obs[:, 4]-correction_median[4,:]) < -27.8)
# print data['snlist'][outlier]

# outlier  = numpy.where((mag_obs[:, 0]-correction_median[0,:]) < -29.35)
# print data['snlist'][outlier]


# print 'Delta', numpy.median(fit['Delta'][:,outlier[0]],axis=0)
# print 'k', numpy.median(fit['k'][:,outlier[0]],axis=0)
# print 'EW0', numpy.median(fit['EW'][:,outlier[0],0],axis=0), EW_renorm[outlier[0],0] 
# print 'EW1', numpy.median(fit['EW'][:,outlier[0],1],axis=0), EW_renorm[outlier[0],1]
fig, axes = plt.subplots()
plt.hist(fit['Delta'].flatten(),normed=True,bins=20)
plt.xlabel(r'$\Delta$')
pp = PdfPages(dirname+'/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['k'].flatten(),normed=True,bins=20)
plt.title(r'$k$')
pp = PdfPages(dirname+'/k_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['k1'].flatten(),normed=True,bins=20)
plt.title(r'$k1$')
pp = PdfPages(dirname+'/k1_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# plt.hist(fit['R'].flatten(),normed=True,bins=20)
# plt.title(r'$D$')
# pp = PdfPages(dirname+'/D_hist.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1]-mag_obs[:,2])
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$B_o-V_o$')
# pp = PdfPages(dirname+'/ebvvsobs.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1],color='red',label=r'$E_\gamma(B-V)$')
# plt.plot(numpy.array([-0.07,0.36]), -29. + 3.96*numpy.array([-0.08,0.36]))
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$B_o$')
# pp = PdfPages(dirname+'/g1.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,0],color='red',label=r'$E_\gamma(B-V)$')
# plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$U_o$')
# pp = PdfPages(dirname+'/g1u.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0), \
#     mag_obs[:,0]- numpy.median(fit['alpha'][:,0][:,None]*fit['EW'][:,:,0]) - numpy.median(fit['beta'][:,0][:,None]*fit['EW'][:,:,1])\
#     ,color='red',label=r'$E_\gamma(B-V)$')
# plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
# plt.xlabel(r'$E_\gamma(B-V)$')
# plt.ylabel(r'$U_o - \alpha_0 EW_{Ca} - \beta_0 EW_{Si}$')
# pp = PdfPages(dirname+'/g2uc.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# plt.scatter(numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],axis=0),mag_obs[:,1]-3.96*numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),color='blue')
# plt.plot(numpy.array([-0.015,0.08]), -29 - 3.46*numpy.array([-0.015,0.08]))
# plt.xlabel(r'$E_\delta(B-V)$')
# plt.ylabel(r'$B_o - \gamma_1 E_\gamma(B-V)$')
# pp = PdfPages(dirname+'/g2.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# plt.hist(((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + 
#   (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1']).flatten(),normed=True,bins=20,
#   label=r'$E(B-V)$')
# plt.xlabel(r'Color Excess')
# # plt.legend()
# pp = PdfPages(dirname+'/greg.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()
# wefe

# plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'], (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']],normed=True,bins=20,
#     label=[r'$E_\gamma(B-V)$',r'$E_{\gamma_1}(B-V)$',r'$E_\delta(B-V)$'],range=(-0.1,0.4))
dustebv = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1']
extebv = (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']
dustebv = dustebv-dustebv[:,0][:,None]
extebv = extebv-extebv[:,0][:,None]
plt.hist([dustebv, extebv],normed=True,bins=20,
  label=[r'$E(B-V)$',r'$E_\delta(B-V)$'])
# plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + 
#   (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']],normed=True,bins=20,
#   label=[r'$E(B-V)$',r'$E_\delta(B-V)$'])
plt.xlabel(r'$E(B-V)-E(B-V)|_0$',fontsize=20)
plt.legend()
pp = PdfPages(dirname+'/ebv.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] ,
  (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1']],normed=True,bins=20,
  label=[r'$E_{\gamma_0}(B-V)$',r'$E_{\gamma_1}(B-V)$'])
plt.xlabel(r'$E(B-V)$')
plt.legend()
pp = PdfPages(dirname+'/ebv_gamma.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist([numpy.median(fit['k'],axis=0) ,numpy.median(fit['k1'],axis=0),numpy.median(fit['R'],axis=0)],normed=True,bins=20,
  label=[r'${\gamma_0}$',r'${\gamma_1}$',r'$R$'])
plt.xlabel(r'$E(B-V)$')
plt.legend()
pp = PdfPages(dirname+'/ebv_natives.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



bins=numpy.arange(-0.3,1,0.05)
crap = fit['gamma'][:,2][:,None]*fit['k']
crap = crap-crap[:,0][:,None]
plt.hist(numpy.median(crap,axis=0),bins,label='median',normed=True,alpha=0.5,width=0.025)
plt.hist(crap.flatten(),bins,label='ideogram',normed=True,alpha=0.5)
plt.xlabel(r'$\gamma^0_2 k_0 - \gamma^0_2 k_0|_0\approx A^F_V|_{R^F=2.44}$')
plt.legend()
plt.tight_layout()
pp = PdfPages(dirname+'/deltagamma0_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

bins=numpy.arange(-0.35,0.1,0.02)
crap2 = fit['gamma1'][:,2][:,None]*fit['k1']
crap2 = crap2-crap2[:,0][:,None]
plt.hist(numpy.median(crap2,axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
plt.hist(crap2.flatten(),bins,label='ideogram',normed=True,alpha=0.5)
plt.xlabel(r'$\gamma^1_2 k_1 - \gamma^1_2 k_1|_0$')
plt.legend(loc=2)
plt.tight_layout()
pp = PdfPages(dirname+'/deltagamma1_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

wefwe

bins = numpy.arange(-0.1,0.121,0.01)
# # plt.hist(fit['Delta'].flatten(),bins,label='ideogram',normed=True,alpha=0.5)
# # plt.hist(numpy.median(fit['Delta'],axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
# # plt.legend()

print 'E_delta(B-V) median sd ',numpy.median(extebv,axis=0).std()
print 'E_ext(B-V) median sd ',numpy.median(dustebv,axis=0).std()

plt.hist( extebv.flatten(),normed=True,bins=bins,alpha=0.5,
  label='ideogram')
plt.hist( numpy.median(extebv,axis=0),normed=True,bins=bins,alpha=0.5,width=0.005,
  label='median')
plt.xlabel(r'$E_\delta(B-V)-E_\delta(B-V)|_0$',fontsize=20)
plt.xlim((-.1,.12))
plt.legend()
pp = PdfPages(dirname+'/ebv_delta.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
# print snnames[numpy.argsort(numpy.median(extebv,axis=0))]
# print numpy.argmin(numpy.median(extebv,axis=0))

plt.hist(numpy.median(extebv[:,1:],axis=0)/numpy.std(extebv[:,1:],axis=0),normed=True,bins=20)
plt.xlabel(r'pull $E_\delta(B-V)$',fontsize=20)
# plt.xlim((-.015,.015))
plt.legend()
pp = PdfPages(dirname+'/ebv_delta_pull.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
# print numpy.where(numpy.median(extebv,axis=0)/numpy.std(extebv,axis=0) < -0.4)
# wefwef
# print snnames[numpy.median(extebv,axis=0)/numpy.std(extebv,axis=0) < -0.4]


# au = fit['gamma'][:,0][:,None]*fit['k'] + fit['rho1'][:,0][:,None]*fit['R']+0.16
# ab = fit['gamma'][:,1][:,None]*fit['k'] + fit['rho1'][:,1][:,None]*fit['R']+0.08
# av = fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R']
# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0)
# (y, ymin, ymax) = numpy.percentile(au-av,(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# x_=numpy.array([-0.03,0.45])
# plt.plot(x_,(4.87-2.96)*x_)
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$E_{eff}(U-V)$')
# pp = PdfPages(dirname+"/Rveff.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (y, ymin, ymax) = numpy.percentile((au-av)/(ab-av),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.ylim((1.8,2.3))
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$R_{eff,U} - R_{eff,V} = E_{eff}(U-V) /E_{eff}(U-V)$')
# pp = PdfPages(dirname+"/Rveff2.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
# (y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$R_{eff, V} = A_{eff, V}/E_{eff}(B-V)$')
# pp = PdfPages(dirname+'/Rveff2.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
# (y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y/(ymax-ymin)*2,xerr=[x-xmin,xmax-x],fmt='o')
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$S/N(R_{eff, V})$')
# pp = PdfPages(dirname+'/Rveff3.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.errorbar(x,y-3.96*x,xerr=[x-xmin,xmax-x],yerr=[(y-ymin),(ymax-y)],fmt='o')
# plt.xlabel(r'$E_{eff}(B-V)+const$')
# plt.ylabel(r'$\Delta A_B  = A_{eff, B}-3.96 (E_{eff}(B-V)+const)$')
# plt.ylim((-0.8,0.1))
# pp = PdfPages(dirname+'/Rveffres.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (y, ymin, ymax) = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] ,(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# (y, ymin, ymax) = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'] ,(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.xlim((-0.1,0.1))
# plt.show()
# (y, ymin, ymax) = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k']- ((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']) ,(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.show()
# x_=numpy.array([-0.08,0.36])
# plt.plot(x_,3.96*x_)
# plt.xlim((-0.1,0.4))
# plt.xlabel(r'$E_{eff}(B-V)+const$')
# plt.ylabel(r'$A_{eff, B}$')
# pp = PdfPages(dirname+'/ebvebv.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()





# IMPORTANT TO KEEP THIS
# import scipy.odr.odrpack as odrpack
# dcolor = numpy.std(ab-av,axis=0)
# dmu = numpy.std(ab,axis=0)

# def f(B, x):
#     return B[0]*x + B[1]
# linear = odrpack.Model(f)
# # mydata = odrpack.Data(x, y, wd=1./np.power(sx,2), we=1./np.power(sy,2))
# mydata = odrpack.RealData(x, y-3.96*x, sx=dcolor, sy=dmu)

# myodr = odrpack.ODR(mydata, linear, beta0=[1., 2.])
# myoutput = myodr.run()
# myoutput.pprint()


# plt.hist(rv)
# pp = PdfPages(dirname+'/Rveff.pdf')
# plt.xlabel(r'$R_{V,eff}$')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median(ab-av,axis=0), numpy.median(rv,axis=0))
# plt.xlabel(r'$E^o(B-V)$')
# plt.ylabel(r'$R_{V,eff}$')
# pp = PdfPages(dirname+'/EVRveff.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
pp = PdfPages(dirname+'/c_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages(dirname+'/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['alpha']-fit['alpha'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$"])
# pp = PdfPages(dirname+'/alpham4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()



figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages(dirname+'/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['eta'],labels=[r"${\eta}_0$",r"${\eta}_1$",r"${\eta}_2$",r"${\eta}_3$",r"${\eta}_4$"])
pp = PdfPages(dirname+'/eta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['beta']-fit['beta'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$"])
# pp = PdfPages(dirname+'/betam4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


figure = corner.corner(fit['gamma'],labels=[r"${\gamma^0}_0$",r"${\gamma^0}_1$",r"${\gamma^0}_2$",r"${\gamma^0}_3$",r"${\gamma^0}_4$"],label_kwargs={'fontsize':22})
pp = PdfPages(dirname+'/gamma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['gamma1'],labels=[r"${\gamma^1}_0$",r"${\gamma^1}_1$",r"${\gamma^1}_2$",r"${\gamma^1}_3$",r"${\gamma^1}_4$"])
pp = PdfPages(dirname+'/gamma1_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = numpy.array([fit['rho11'],fit['rho12'],fit['rho13']])
# mega = numpy.transpose(mega)
# figure = corner.corner(mega,labels=[r"${\delta}_{0}$",r"${\delta}_{1}$",r"${\delta}_{2}$"])
# pp = PdfPages(dirname+'/rho_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


figure = corner.corner(fit['rho1'],labels=[r"${\delta}_{U}$",r"${\delta}_{B}$",r"${\delta}_{V}$",r"${\delta}_{R}$",r"${\delta}_{I}$"],\
  truths=[0,0,0,0,0])
pp = PdfPages(dirname+'/delta_corner.pdf')
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# wmn, wless, wmore = numpy.percentile(fit['rho1'][:,0][:,None]*fit['R'],(50,50-34,50+34),axis=0)
# shit = wmn/(wmn-wless)
# shit[wmn < 0] = wmn[wmn < 0]/(wmn[wmn < 0]-wmore[wmn < 0])
# wmax = numpy.argmax(shit)
# temprho1 = numpy.array(fit['rho1'])
# for i in xrange(temprho1.shape[0]):
#     # if (scipy.stats.skew(fit['R'][i,:]) < 0):
#     if (fit['R'][i,wmax] < 0):
#         temprho1[i,:] *=-1

temprho1,tempD = flip.flip(fit)

m1,m2 = numpy.percentile(fit['R']*fit['rho1'][:,0][:,None],(50-34,50+34),axis=0)
print numpy.median((m2-m1)/2)

fig, axes = plt.subplots()
plt.hist(numpy.median(fit['R']*fit['rho1'][:,0][:,None],axis=0),bins=20)
plt.xlabel(r'$A_{\delta U}$')
pp = PdfPages(dirname+'/AU_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(temprho1,labels=[r"${\delta}_{U}$",r"${\delta}_{B}$",r"${\delta}_{V}$",r"${\delta}_{R}$",r"${\delta}_{I}$"],\
  truths=[0,0,0,0,0])
pp = PdfPages(dirname+'/delta_corner_flipped.pdf')
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# Rsd = numpy.linalg.norm(fit['rho1'],axis=1)
# figure = corner.corner(fit['rho1']/Rsd[:,None],labels=[r"${\delta}_{U}$",r"${\delta}_{B}$",r"${\delta}_{V}$",r"${\delta}_{R}$",r"${\delta}_{I}$"],\
#   truths=[0,0,0,0,0])
# pp = PdfPages(dirname+'/delta_corner_scaled.pdf')
# plt.tight_layout()
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# figure = corner.corner(fit['rho1']/numpy.linalg.norm(fit['rho1'],axis=1)[:,None],labels=[r"${\delta}_{U}$",r"${\delta}_{B}$",r"${\delta}_{V}$",r"${\delta}_{R}$",r"${\delta}_{I}$"],\
#   truths=[0,0,0,0,0])
# pp = PdfPages(dirname+'/delta_unit_corner.pdf')
# plt.tight_layout()
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

dum = (fit['rho1']-fit['rho1'][:,2][:,None])/((fit['rho1'][:,0]-fit['rho1'][:,1])[:,None])
dum = numpy.delete(dum,2,axis=1)
figure = corner.corner(dum,range=((0,2),(-1,1),(-0.5,1.5),(0,2)))
pp = PdfPages(dirname+'/g_corner.pdf')
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



figure = corner.corner(fit['rho1'][:,1:]/fit['rho1'][:,0][:,None],labels=[r"${\delta}_{0}/{\delta}_{U}$",r"${\delta}_{1}/{\delta}_{U}$",\
  r"${\delta}_{2}/{\delta}_{U}$",r"${\delta}_{3}/{\delta}_{U}$"], \
  range=[[-3.5,3.5],[-3.5,3.5],[-3.5,3.5],[-3.5,3.5]])
pp = PdfPages(dirname+'/deltaratio_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None]

# mega = mega-mega[:,2][:,None]
# mega = numpy.delete(mega,2,axis=1)
# mega = numpy.delete(mega,1,axis=1)
# figure = corner.corner(mega,labels=[r"$R_0$",r"$R_3$",r"$R_4$"])
# pp = PdfPages(dirname+'/rx-rv_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages(dirname+'/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# plt.hist(fit['lp__'],cumulative=True) #
# plt.xlabel(r'$\ln{p}$')
# pp = PdfPages(dirname+'/lp.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


with PdfPages(dirname+'/multipage_pdf.pdf') as pdf:

    lineobjects = plt.plot(fit['lp__'][::10])
    plt.title(r'log p')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(fit['L_sigma'][::10,:],label=[])
    plt.title(r'$L_\sigma$')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
    
    plt.plot(fit['EW'][::10,:10,0])
    plt.title('EW - 0')
    pdf.savefig()
    plt.close()

    plt.plot(fit['EW'][::10,:10,1])
    plt.title('EW - 1')
    pdf.savefig()
    plt.close()



    plt.plot(fit['Delta_scale'][::10])
    plt.title(r'$\Delta$ scale')
    pdf.savefig()
    plt.close()
    
    lineobjects = plt.plot(fit['alpha'][::10],alpha=0.5,label=['U','B','V','R','I'])
    plt.title('alpha')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
    
    lineobjects = plt.plot(fit['beta'][::10],alpha=0.5,label=['U','B','V','R','I'])
    plt.title('beta')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
    
    lineobjects = plt.plot(fit['gamma'][::10],label=['U','B','V','R','I'])
    plt.title('gamma')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
 
    mega = fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None]

    figure = corner.corner(mega,labels=[r"$R_0$",r"$R_1$",r"$R_2$",r"$R_3$",r"$R_4$"])
    pp = PdfPages(dirname+'/rx_corner.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()

    lineobjects = plt.plot(fit['gamma1'][::10],label=['U','B','V','R','I'])
    plt.title('gamma1')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
 
    # mega = fit['gamma1']/((fit['gamma1'][:,1]-fit['gamma1'][:,2]))[:,None]

    # figure = corner.corner(mega,labels=[r"$R_0$",r"$R_1$",r"$R_2$",r"$R_3$",r"$R_4$"])
    # pp = PdfPages(dirname+'/r1x_corner.pdf')
    # plt.savefig(pp,format='pdf')
    # pp.close()
    # plt.close()

    # mega = fit['rho1']/((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None]

    # figure = corner.corner(mega,labels=[r"$R_{\delta 0}$",r"$R_{\delta 1}$",r"$R_{\delta 2}$",r"$R_{\delta 3}$",r"$R_{\delta 4}$"], \
    #     range=[[-8.,0.5] for x in xrange(5)])
    # pp = PdfPages(dirname+'/rxdelta_corner.pdf')
    # plt.savefig(pp,format='pdf')
    # pp.close()
    # plt.close()

    # plt.plot(fit['R'][::10,:10])
    # plt.title('R')
    # pdf.savefig()
    # plt.close()

    # mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['L_sigma']])
    # mega = numpy.transpose(mega)

    # for index in xrange(5):
    #     figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index),r"$\alpha_{}$".format(index),\
    #         r"$\beta_{}$".format(index),r"$\sigma_{}$".format(index)])
    #     pdf.savefig()
    #     plt.close()

    # filter wavelengths
# r = numpy.log(8635./3276.)
# edges= 3276.*numpy.exp(numpy.arange(6)*r/5.)
# # edges = numpy.array([[3300., 4102], [4102, 5100], [5200, 6289], [6289, 7607], [7607, 9200]])
# efflam = []
# for i in xrange(5):
#     efflam.append((edges[i+1]-edges[i])/2+edges[i])
#     # efflam.append((edges[i][1]-edges[i][0])/2+edges[i][0])

#from manu
efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])
# [3701, 4601, 5744, 6948, 8403]
rc('text', usetex=True)



filts = ['U','B','V','R','I']
# lambdas = numpy.arange(3000.,9000,100)
# rvs=[2.97]
# avs = [12.25]
# for rv,av in zip(rvs,avs):
#   f99 = sncosmo.F99Dust(r_v =rv)
#   f99.set(ebv=av)
#   A_ = f99.propagate(lambdas,1)
#   A_=-2.5*numpy.log10(A_)

#   # norm = f99.propagate([efflam[1]],1.) / f99.propagate(numpy.array([efflam[2]]),1.)
#   # norm  = -2.5*numpy.log10(norm)
#   # A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
#   # norm  = sncosmo._extinction.ccm89(numpy.array([efflam[1]]), 1., rv)-sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
#   # A_ = A_/norm[0]
#   plt.plot(lambdas,A_,label=r"$R^F_V={:.2f}$".format(rv))

# (y, ymin, ymax) = numpy.percentile(fit['gamma'],(50,50-34,50+34),axis=0)
# plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
# plt.legend()
# plt.xlabel(r'Wavelength (\AA)')
# plt.ylabel(r'$\gamma_X$')
# pp = PdfPages(dirname+'/fitz.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


nlinks = fit['gamma'].shape[0]
mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['gamma1'],fit['rho1'],fit['L_sigma']])
mega = numpy.transpose(mega)


cname=['U','B','V','R','I']

for index in xrange(5):
    figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(cname[index]), r"$\alpha_{}$".format(cname[index]),\
                    r"$\beta_{}$".format(cname[index]),r"$\eta_{}$".format(cname[index]),r"$\gamma^0_{}$".format(cname[index]),\
                    r"$\gamma^1_{}$".format(cname[index]),r"$\delta_{{{}}}$".format(cname[index]), r"$\sigma_{}$".format(cname[index])], \
                    truths=[None,0,0,0,0,0,0,0],label_kwargs={'fontsize':22})
    figure.suptitle(filts[index],fontsize=28)

    for ax in figure.get_axes():
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

    pp = PdfPages(dirname+'/coeff{}.pdf'.format(index))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()



lambdas = numpy.arange(3000.,9000,100)
# for rv in rvs:

#     f99 = sncosmo.F99Dust(r_v =rv)
#     f99.set(ebv=1.37)
#     A_ = f99.propagate(lambdas,1.)
#     A_=-2.5*numpy.log10(A_)

#     norm = f99.propagate([efflam[1]],1.) / f99.propagate(numpy.array([efflam[2]]),1.)
#     norm  = -2.5*numpy.log10(norm)
#     # A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
#     # norm  = sncosmo._extinction.ccm89(numpy.array([efflam[1]]), 1., rv)-sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
#     A_ = A_/norm[0]
#     plt.plot(lambdas,A_,label=r"$R^F_V={:.2f}$".format(rv))

# (y, ymin, ymax) = numpy.percentile(fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]),(50,50-34,50+34),axis=0)
# plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
# plt.legend()
# plt.xlabel(r'Wavelength (\AA)')
# plt.ylabel(r'$R_X=\frac{\gamma_X}{\gamma_1-\gamma_2}$')
# pp = PdfPages(dirname+'/ccm.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# for rv in rvs:

#     f99 = sncosmo.F99Dust(r_v =rv)
#     f99.set(ebv=1.37)
#     A_ = f99.propagate(lambdas,1.)

#     norm  = f99.propagate(numpy.array([efflam[2]]),1.) 
#     # A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
#     # norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
#     # A_ = A_/norm[0]
#     A_ = numpy.log10(A_)/numpy.log10(norm)
#     plt.plot(lambdas,A_,label=r"$R^F_V={:.1f}$".format(rv))

# (y, ymin, ymax) = numpy.percentile(fit['gamma']/fit['gamma'][:,2][:,None],(50,50-34,50+34),axis=0)
# plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
# plt.legend()
# plt.xlabel(r'Wavelength (\AA)')
# plt.ylabel(r'$\frac{\gamma_X}{\gamma_2}$')
# pp = PdfPages(dirname+'/ccm2.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()



labels = list('UBVRI')
from matplotlib.ticker import FuncFormatter, MaxNLocator
def format_fn2(tick_val, tick_pos):
    if int(tick_val) in numpy.arange(5):
        return labels[int(tick_val)]
    else:
        return ''
(y, ymin, ymax) = numpy.percentile(fit['rho1']/fit['rho1'][:,0][:,None],(50,50-34,50+34),axis=0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(numpy.arange(5),y,yerr=[y-ymin,ymax-y],fmt='o')
ax.xaxis.set_major_formatter(FuncFormatter(format_fn2))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.axhline(0,linestyle=':')
ax.set_xlim((-0.5,4.5))
ax.set_xlabel(r'Band $X$')
ax.set_ylabel(r'$\frac{\delta_X}{\delta_U}$')
ax.set_ylim((-1.5,1.2))
plt.tight_layout()
pp = PdfPages(dirname+'/deltaratio.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(y, ymin, ymax) = numpy.percentile(temprho1/temprho1[:,0][:,None],(50,50-34,50+34),axis=0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(numpy.arange(5),y,yerr=[y-ymin,ymax-y],fmt='o')
ax.xaxis.set_major_formatter(FuncFormatter(format_fn2))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.axhline(0,linestyle=':')
ax.set_xlim((-0.5,4.5))
ax.set_xlabel(r'Band $X$')
ax.set_ylabel(r'$\frac{\delta_X}{\delta_U}$')
ax.set_ylim((-1.3,1.3))
plt.tight_layout()
pp = PdfPages(dirname+'/deltaratio_flip.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), \
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']).flatten(),\
    ((fit['gamma1'][:,1] - fit['gamma1'][:,2])[:,None]*fit['k1']).flatten(),\
    ((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']).flatten()])

mega = numpy.transpose(mega)
mega=mega[::50,:]

figure = corner.corner(mega,labels=[r"$\Delta$",r"$EW_{Ca}$",r"$EW_{Si}$",r"$\lambda_{Si}$",r"$E_{\gamma^0}(B-V)$",r"$E_{\gamma^1}(B-V)$",r"$E_\delta(B-V)$"], \
    range=numpy.zeros(7)+1.,label_kwargs={'fontsize':22})
for ax in figure.get_axes():
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
pp = PdfPages(dirname+'/perobject_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()





