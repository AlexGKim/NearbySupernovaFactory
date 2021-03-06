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
from chainconsumer import ChainConsumer

mpl.rcParams['font.size'] = 28
# rc('text', usetex=True)

f = open('fix3_x1_decorr.pkl','rb')
# (fit,_) = pickle.load(f)
fit = pickle.load(f)

# for key in fit.keys():
#     print key, fit[key].min(), fit[key].max()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

# dic_phreno=cPickle.load(open("phrenology_2016_12_01_CABALLOv1.pkl"))

# dic_meta=cPickle.load(open("META.pkl"))


sivel,sivel_err,x1,x1_err, zcmb, _, _, _, c = sivel.sivel(data)
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
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]
snname=numpy.array(data['snlist'])[use]

# dum = [[snnm,z_] for snnm,z_ in zip(snname, zcmb) if (z_ < 0.025 or z_ > 0.085)]
# dum = numpy.array(dum)
# inde = numpy.argsort(dum[:,1])
# for i in inde:
#     print dum[i,0], dum[i,1]



dic_meta=cPickle.load(open("META.pkl"))
dic_meta['PTF12efn']={'host.zhelio':0.06194}
# dic_meta['PTF12efn']['host.zhelio']=0.06194
zhelio=[]
for nm in snname:
    zhelio.append(dic_meta[nm]['host.zhelio'])
zhelio=numpy.array(zhelio)

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


trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))


# correction = [fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
#     + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']\
#     for i in xrange(5)]
# correction = numpy.array(correction)

# intrinsic = numpy.median(correction,axis=1)
# fixev = fit['ev']
# fixev = fixev * numpy.sign(fixev[:,2])[:,None]
# corner.corner(fit['ev_sig'][:,None]*fixev,labels=[r"$\sigma_p \phi_{\hat{U}}$",r"$\sigma_p \phi_{\hat{B}}$",r"$\sigma_p \phi_{\hat{V}}$",r"$\sigma_p \phi_{\hat{R}}$",r"$\sigma_p \phi_{\hat{I}}$"])
# pp = PdfPages('output_fix3_x1/sigev.pdf')
# plt.savefig(pp,format='pdf',bbox_inches='tight')
# pp.close()
# plt.close()


#significance of ev
plt.hist(fit['ev_sig'])
plt.xlabel(r'$\sigma_p$')
pp = PdfPages('output_fix3_x1/sigevhist.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()

# output = open('intrinsic.pkl','wb')
# print intrinsic.shape, numpy.array(data['snlist'])[use].shape
# pickle.dump([intrinsic, numpy.array(data['snlist'])[use]], output, protocol=2)
# output.close()
# wefew
# colors = fit['c']+mag_mn[None,:]
# colors  = colors - colors[:,2][:,None]
# colors = numpy.delete(colors, 2,axis=1)
# figure = corner.corner(colors,labels=[r"$U_0-V_0$",r"$B_0-V_0$",r"$R_0-V_0$",r"$I_0-V_0$"])
# pp = PdfPages('output_fix3_x1/col_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# plt.hist(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),bins=20)
# plt.xlabel(r'$E_{\gamma_0}(B-V)$')
# plt.legend()
# pp = PdfPages('output_fix3_x1/ebv_gamma0.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# bins = numpy.arange(-0.2,0.8001,0.02)
# plt.hist(fit['Delta'].flatten(),bins,label='ideogram',normed=True,alpha=0.5)
# plt.hist(numpy.median(fit['Delta'],axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)

bins=numpy.arange(-0.3,1,0.05)
plt.hist(numpy.median((fit['gamma'][:,2])[:,None]*fit['k'],axis=0),bins,label='median',normed=True,alpha=0.5,width=0.02)
plt.hist(((fit['gamma'][:,2])[:,None]*fit['k']).flatten(),bins,label='ideogram',normed=True,alpha=0.5)
plt.xlabel(r'$\gamma^0_2 k_0 \approx A^F_V|_{R^F=2.44}$')
plt.legend()

plt.tight_layout()
pp = PdfPages('output_fix3_x1/gamma0_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.clf()
bins=numpy.arange(-0.2,1.2,0.05)
crap = fit['gamma'][:,2][:,None]*fit['k']
crap = crap-crap[:,0][:,None]
plt.hist(crap.flatten(),bins,label='posterior stack',normed=True,alpha=0.5)
plt.hist(numpy.median(crap,axis=0),bins,label='median',normed=True,alpha=0.5,width=0.025)
plt.xlabel(r'$\gamma^0_{\hat{V}} g_0 - \gamma^0_{\hat{V}} g_0|_0 \approx A^F_V|_{\langle R^F_{eff}\rangle=2.43}$')
plt.legend(fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
pp = PdfPages('output_fix3_x1/deltagamma0_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.clf()
bins=numpy.arange(-0.4,0.2,0.02)
crap2 = fit['rho1'][:,2][:,None]*fit['R']
crap2 = crap2-crap2[:,0][:,None]
plt.hist(crap2.flatten(),bins,label='posterior stack',normed=True,alpha=0.5)
plt.hist(numpy.median(crap2,axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
plt.xlabel(r'$\gamma^1_{\hat{V}} g_1 - \gamma^1_{\hat{V}} g_1|_0$')
plt.legend(loc=2,fontsize=20)
plt.xlim((-0.3,0.1))
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
pp = PdfPages('output_fix3_x1/deltagamma1_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# #conversion to Fitzpatrick parameters A and E(B-V)
# tmat= numpy.array([[2.82,1.15],[-5.27,0.72]])
# tmat = tmat.T
# shit=[]
# for i in xrange(fit['k'].shape[1]):
#     knorm =  numpy.median(fit['k'][:,i]*(fit['gamma'][:,1]-fit['gamma'][:,2]))
#     Rnorm = numpy.median(fit['R'][:,i]*(fit['rho1'][:,1]-fit['rho1'][:,2]))
#     shit.append(numpy.dot(tmat,numpy.array([knorm,Rnorm])))
# shit=numpy.array(shit)
# print shit.shape
# plt.scatter(shit[:,1],shit[:,0])
# plt.plot([-.1,.5],[-2.44*.1,2.44*.5])
# plt.show()

# wefwe

bins=numpy.arange(-0.35,0.1,0.02)
plt.hist(numpy.median((fit['rho1'][:,2])[:,None]*fit['R'],axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
plt.hist(((fit['rho1'][:,2])[:,None]*fit['R']).flatten(),bins,label='ideogram',normed=True,alpha=0.5)
plt.xlabel(r'$\gamma^1_2 k_1$')
plt.legend(loc=2)
plt.tight_layout()
pp = PdfPages('output_fix3_x1/gamma1_med.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(x, xmin, xmax) = numpy.percentile(crap,(50,50-34,50+34),axis=0)
(y, ymin, ymax) = numpy.percentile(crap2,(50,50-34,50+34),axis=0)
plt.errorbar(x, y, xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-ymin],fmt='o',alpha=0.2)
plt.scatter(x, y, s=1,alpha=0.8)
plt.xlabel(r'$\gamma^0_V k_0 - \gamma^0_V k_0|_0\approx A^F_V|_{R^F=2.44}$')
plt.ylabel(r'$\gamma^1_V k_1 - \gamma^1_V k_1|_0$')
pp = PdfPages("output_fix3_x1/deltakk.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(x, xmin, xmax) = numpy.percentile(fit['gamma'][:,2][:,None]*fit['k'],(50,50-34,50+34),axis=0)
(y, ymin, ymax) = numpy.percentile(fit['rho1'][:,2][:,None]*fit['R'],(50,50-34,50+34),axis=0)
plt.errorbar(x, y, xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-ymin],fmt='o',alpha=0.2)
plt.scatter(x, y, s=1,alpha=0.8)
plt.xlabel(r'$\gamma^0_V k_0 \approx A^F_V|_{R^F=2.44}$')
plt.ylabel(r'$\gamma^1_V k_1$')
pp = PdfPages("output_fix3_x1/kk.pdf")
plt.tight_layout()
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# crap=fit['rho1'][:,2][:,None]*fit['R']
# med = numpy.percentile(x,(25,75))
# plt.hist([crap[:,x<=med[0]].flatten(),crap[:,x>med[1]].flatten()],label=[r'Low $\gamma^0_V k_0$', r'High $\gamma^0_V k_0$'])
# plt.legend(loc=2)
# # plt.show()
# plt.clf()

# kappa1=1/2.4
# kappa2 = -6.8
# av = ((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k']+ \
#     kappa2*(fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'])
# e = (kappa1*(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] +\
#     (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'])
# r = av/e

# (x, xmin, xmax) = numpy.percentile(av,(50,50-34,50+34),axis=0)
# (y, ymin, ymax) = numpy.percentile(e,(50,50-34,50+34),axis=0)
# plt.errorbar(x, y, xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-ymin],fmt='o')
# plt.xlabel(r'$\gamma^0_2 k_0$')
# plt.ylabel(r'$\gamma^1_2 k_1$')
# pp = PdfPages("output_fix3_x1/ra.pdf")
# plt.tight_layout()
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# wefwef

plt.hist(numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],axis=0),bins=20)
plt.xlabel(r'$E_{\gamma_1}(B-V)$')
plt.legend()
pp = PdfPages('output_fix3_x1/ebv_gamma1.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(numpy.concatenate((EW_obs,sivel[:,None]),axis=1),labels=[r"$EW_{Ca,o}$",r"$EW_{Si,o}$",r"$\lambda_{Si,o}$"])
for ax in figure.get_axes():
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        
pp = PdfPages('output_fix3_x1/feature_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(x, xmin, xmax) = numpy.percentile( ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']),(50,50-34,50+34),axis=0)
(y, ymin, ymax) = numpy.percentile(((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']),(50,50-34,50+34),axis=0)

plt.errorbar(x,y,xerr=(x-xmin,xmax-x), yerr=(y-ymin,ymax-y),fmt='o')
plt.xlim((-0.08,0.05))
# plt.ylim((-0.02,0.05))
plt.xlabel(r'$E_{\gamma^0}(B-V)$')
plt.ylabel(r'$E_{\gamma^1}(B-V)$')
pp = PdfPages("output_fix3_x1/egammaedelta_corner.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


(y, ymin, ymax) = numpy.percentile(fit['EW'][:,:,1],(50,50-34,50+34),axis=0)
plt.errorbar(x1, y, xerr=[x1_err,x1_err],yerr=[y-ymin,ymax-ymin],fmt='o')
plt.xlabel(r'$X_1$')
plt.ylabel(r'$EW_{Si}$')
pp = PdfPages("output_fix3_x1/x1si.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(y, ymin, ymax) = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'][:,:],(50,50-34,50+34),axis=0)
plt.errorbar(x1, y, xerr=[x1_err,x1_err],yerr=[y-ymin,ymax-ymin],fmt='o')
plt.xlabel(r'$X_1$')
plt.ylabel(r'$E_\delta(B-V)$')
pp = PdfPages("output_fix3_x1/x1D.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


offset = 0.01
au = fit['gamma'][:,0][:,None]*fit['k'] + fit['rho1'][:,0][:,None]*fit['R']
ab = fit['gamma'][:,1][:,None]*fit['k'] + fit['rho1'][:,1][:,None]*fit['R']
av = fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R']


plt.hist(numpy.median(av,axis=0),bins=20)
plt.xlabel(r'$A_V$')
pp = PdfPages("output_fix3_x1/Av_med.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0)
(y, ymin, ymax) = numpy.percentile(au-av,(50,50-34,50+34),axis=0)

plt.scatter(x ,numpy.median(av/(ab-av),axis=0))
plt.ylim((-2,6))
plt.xlim((-0.1,.4))
plt.ylabel(r'$R_{eff,V}$')
plt.xlabel(r'$E_{eff}(B-V)$')
pp = PdfPages("output_fix3_x1/Rveff4.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


x=x-x.min()-offset*1.3
y=y-y.min()-offset
plt.scatter(x,y)
# x_=numpy.array([-0.03,0.45])
# plt.plot(x_,(4.87-2.96)*x_)
plt.xlim((-0.04,0.5))
plt.ylim((-0.04,1))
plt.xlabel(r'$E_{eff}(B-V)$')
plt.ylabel(r'$E_{eff}(U-V)$')
pp = PdfPages("output_fix3_x1/Rveff.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(y, ymin, ymax) = numpy.percentile((au-av-y.min()-offset)/(ab-av-x.min()-offset),(50,50-34,50+34),axis=0)
plt.scatter(x,y)
# plt.ylim((1.8,2.3))
#plt.xlim((0,0.5))
plt.xlabel(r'$E_{eff}(B-V)$')
plt.ylabel(r'$R_{eff,U} - R_{eff,V}$')
pp = PdfPages("output_fix3_x1/Rveff2.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

temp = (fit['gamma'][:,1][:,None]*fit['k']-fit['gamma'][:,2][:,None]*fit['k'] ) /\
 (fit['rho1'][:,1][:,None]*fit['R'] - fit['rho1'][:,2][:,None]*fit['R'])
(x, xmin, xmax) = numpy.percentile(temp,(50,50-34,50+34),axis=0)
plt.scatter(x,y)
# plt.ylim((1.8,2.3))
plt.xlabel(r'$E_\gamma(B-V)/E_\rho(B-V)$')
plt.ylabel(r'$R_{eff,U} - R_{eff,V}$')
pp = PdfPages("output_fix3_x1/Rveff3.pdf")
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


filts = ['U','B','V','R','I']
correction = [fit['Delta']+ fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
    + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']+ fit['zeta'][:,i][:,None]*fit['x1']\
    + fit['gamma'][:,i][:,None]*fit['k']+ fit['rho1'][:,i][:,None]*fit['R'] + (fit['ev_sig']*fit['ev'][:,i])[:,None]* fit['mag_int_raw'] \
    for i in xrange(5)]

correction = numpy.array(correction)

# correction subtracting out the zero
correction_dot = correction-correction[:,:,0][:,:,None]
correction_dot = correction_dot[:,:,1:]

correction_dot_mn = correction_dot.mean(axis=1)
correction_dot_std = correction_dot.std(axis=1)
correction_dot_mn = numpy.swapaxes(correction_dot_mn,0,1)
correction_dot_std = numpy.swapaxes(correction_dot_std,0,1)

colcorrection = numpy.zeros((correction.shape[1],correction.shape[2],5,5))

for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        colcorrection[:,:,m0,m1] = correction[m0,:,:] - correction[m1,:,:]

colcorrection_dot = colcorrection[:,:,:,:] - colcorrection[:,0,:,:][:,None,:,:]
colcorrection_dot = colcorrection_dot[:,1:,:,:]
colcorrection_dot_mn = colcorrection_dot.mean(axis=0)
colcorrection_dot_std = colcorrection_dot.std(axis=0)

mag_renorm_dot = mag_renorm[:,:] -mag_renorm[0,:]
mag_renorm_dot  = mag_renorm_dot[1:,:]

mag_renorm_dot_var = [numpy.diag(mag_cov[index]) for index in xrange(nsne)]
mag_renorm_dot_var = mag_renorm_dot_var + mag_renorm_dot_var[0]
mag_renorm_dot_var = mag_renorm_dot_var[1:]

obscolors = numpy.zeros((nsne,5,5))
obscolors_var = numpy.zeros((nsne,5,5))
for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        obscolors[:,m0,m1]=mag_renorm[:,m0]-mag_renorm[:,m1]
        obscolors_var[:,m0,m1]=mag_cov[:,m0,m0]+ mag_cov[:,m1,m1]-2*mag_cov[:,m0,m1]
obscolors_dot = obscolors[:,:,:] - obscolors[0,:,:]
obscolors_dot = obscolors[1:,:,:]
obscolors_dot_var = obscolors_var[:,:,:] + obscolors_var[0,:,:]
obscolors_dot_var = obscolors_dot_var[1:,:,:]

pp = PdfPages("mpull.pdf")
magresidual  = mag_renorm_dot - correction_dot_mn
magresidual_std = numpy.sqrt(correction_dot_std**2 + mag_renorm_dot_var)
magresidual_mn = numpy.sum(magresidual/magresidual_std**2,axis=0)/ numpy.sum(1/magresidual_std**2,axis=0)

magresidual_pull = (magresidual- magresidual_mn[None,:])/magresidual_std

colresidual = obscolors_dot - colcorrection_dot_mn
colresidual_std = numpy.sqrt(colcorrection_dot_std**2 + obscolors_dot_var)
colresidual_mn = numpy.sum(colresidual/colresidual_std**2,axis=0)/ numpy.sum(1/colresidual_std**2,axis=0)
colresidual_pull = (colresidual- colresidual_mn[None,:])/colresidual_std

for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        fig, axes = plt.subplots(nrows=5,sharex=True)
        for i in xrange(5):
            axes[i].errorbar(obscolors_dot[:,m0,m1],magresidual[:,i], \
                xerr=[numpy.sqrt(obscolors_dot_var[:,m0,m1]),numpy.sqrt(obscolors_dot_var[:,m0,m1])], \
                yerr=[magresidual_std[:,i],magresidual_std[:,i]],\
                linestyle='None',alpha=0.5,fmt='o')
            axes[i].set_ylabel(r"$\hat{{{}}}$".format(filts[i]),fontsize=12)
            for tick in axes[i].yaxis.get_major_ticks():
                    tick.label.set_fontsize(8)
            axes[i].set_ylim((-0.2,0.2))
            # axes[i].text(0.05, 0.9, \
            #     "RMS: Pull {:5.2f} Mag {:6.3f}".format(numpy.std(magresidual_pull[:,i]),numpy.std(magresidual[:,i])), \
            #     horizontalalignment='left', verticalalignment='center', \
            #     transform=axes[i].transAxes,fontsize=8)           
        fig.set_size_inches(8,11)
        axes[4].set_xlabel(r"$\hat{{{}}}_o-\hat{{{}}}_o$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[4].xaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
        axes[m0].set_axis_bgcolor('0.8')
        axes[m1].set_axis_bgcolor('0.8')
        plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()

pp = PdfPages("cpull.pdf")

for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        fig, axes = plt.subplots(nrows=5, ncols=2,sharex=True)
        i = 0
        for m0_ in xrange(5):
            for m1_ in xrange(m0_+1,5):
                indeces  = numpy.unravel_index(i,axes.shape)
                axes[indeces[0],indeces[1]].errorbar(obscolors_dot[:,m0,m1], \
                    colresidual[:,m0_,m1_], \
                    xerr=[numpy.sqrt(obscolors_dot_var[:,m0,m1]),numpy.sqrt(obscolors_dot_var[:,m0,m1])], \
                    yerr=[colresidual_std[:,m0_,m1_],colresidual_std[:,m0_,m1_]], \
                    linestyle='None',alpha=0.5,fmt='o')
                axes[indeces[0],indeces[1]].set_ylabel(r"$\hat{{{}}}-\hat{{{}}}$ residual".format(filts[m0_],filts[m1_]),fontsize=12)
                for tick in axes[indeces[0],indeces[1]].yaxis.get_major_ticks():
                        tick.label.set_fontsize(8)
                if m0_==m0 or m0_==m1 or m1_==m0 or m1_==m1:
                    axes[indeces[0],indeces[1]].set_axis_bgcolor('0.8')
                # axes[indeces[0],indeces[1]].text(0.05, 0.9, \
                #     "RMS: Pull {:5.2f} Color {:6.3f}".format(numpy.std(colresidual_pull[:,m0,m1]),numpy.std(colresidual[:,m0,m1])), horizontalalignment='left', verticalalignment='center', \
                #     transform=axes[indeces[0],indeces[1]].transAxes,fontsize=8)
                i=i+1        
        fig.set_size_inches(8,11)
        axes[4,0].set_xlabel(r"$\hat{{{}}}_o-\hat{{{}}}_o$".format(filts[m0],filts[m1]),fontsize=12)
        axes[4,1].set_xlabel(r"$\hat{{{}}}_o-\hat{{{}}}_o$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[4,0].xaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
        for tick in axes[4,1].xaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
        fig.subplots_adjust(hspace=.1, wspace=.27)
        plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


pp = PdfPages("extrinsic.pdf")
fig, axes = plt.subplots(nrows=5, ncols=2)

dust= [fit['gamma'][:,i][:,None]*fit['k']+ fit['rho1'][:,i][:,None]*fit['R'] for i in xrange(5)]
dust=numpy.array(dust)

dustcorrection  = numpy.zeros((dust.shape[1],dust.shape[2],5,5))

for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        dustcorrection[:,:,m0,m1] = dust[m0,:,:] - dust[m1,:,:]

dustcorrection_dot = dustcorrection[:,:,:,:] - dustcorrection[:,0,:,:][:,None,:,:]
dustcorrection_dot = dustcorrection_dot[:,1:,:,:]
dustcorrection_dot_mn = dustcorrection_dot.mean(axis=0)
dustcorrection_dot_std = dustcorrection_dot.std(axis=0)
i=0
for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        indeces  = numpy.unravel_index(i,axes.shape)

        # pfit, V  = numpy.polyfit(numpy.mean(dust[m0,:,:]-dust[m1,:,:],axis=0), residual/colorerror_,1,cov=True)
        axes[indeces[0],indeces[1]].errorbar(dustcorrection_dot_mn[:,m0,m1],colresidual[:,m0,m1], \
            yerr=[colresidual_std[:,m0,m1],colresidual_std[:,m0,m1]],\
            xerr=[dustcorrection_dot_std[:,m0,m1],dustcorrection_dot_std[:,m0,m1]],linestyle='None',alpha=0.5,fmt='o')

        # axes[indeces[0],indeces[1]].plot(axes[indeces[0],indeces[1]].get_xlim(),numpy.array(axes[indeces[0],indeces[1]].get_xlim())*pfit[0]+pfit[1])
        axes[indeces[0],indeces[1]].set_ylabel(r"$\hat{{{}}}-\hat{{{}}}$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[indeces[0],indeces[1]].yaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
            # axes[i].set_ylim((-0.15,0.15))            
        fig.set_size_inches(8,11)
        axes[indeces[0],indeces[1]].set_xlabel(r"Extrinsic $\hat{{{}}}-\hat{{{}}}$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[indeces[0],indeces[1]].xaxis.get_major_ticks():
                tick.label.set_fontsize(8)

        i=i+1
fig.subplots_adjust(hspace=.4, wspace=.22)
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()

pp = PdfPages("intrinsic.pdf")
fig, axes = plt.subplots(nrows=5, ncols=2)



intrinsic= [(fit['ev_sig']*fit['ev'][:,i])[:,None]* fit['mag_int_raw'] for i in xrange(5)]
intrinsic=numpy.array(intrinsic)

intrinsiccorrection  = numpy.zeros((intrinsic.shape[1],intrinsic.shape[2],5,5))

for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        intrinsiccorrection[:,:,m0,m1] = intrinsic[m0,:,:] - intrinsic[m1,:,:]

intrinsiccorrection_dot = intrinsiccorrection[:,:,:,:] - intrinsiccorrection[:,0,:,:][:,None,:,:]
intrinsiccorrection_dot = intrinsiccorrection_dot[:,1:,:,:]
intrinsiccorrection_dot_mn = intrinsiccorrection_dot.mean(axis=0)
intrinsiccorrection_dot_std = intrinsiccorrection_dot.std(axis=0)


i=0
for m0 in xrange(5):
    for m1 in xrange(m0+1,5):
        indeces  = numpy.unravel_index(i,axes.shape)
        axes[indeces[0],indeces[1]].errorbar(intrinsiccorrection_dot_mn[:,m0,m1],colresidual[:,m0,m1], \
            yerr=[colresidual_std[:,m0,m1],colresidual_std[:,m0,m1]],\
            xerr=[intrinsiccorrection_dot_std[:,m0,m1],intrinsiccorrection_dot_std[:,m0,m1]],linestyle='None',alpha=0.5,fmt='o')
        axes[indeces[0],indeces[1]].set_ylabel(r"$\hat{{{}}}-\hat{{{}}}$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[indeces[0],indeces[1]].yaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
            # axes[i].set_ylim((-0.15,0.15))            
        fig.set_size_inches(8,11)
        axes[indeces[0],indeces[1]].set_xlabel(r"Intrinsic $\hat{{{}}}-\hat{{{}}}$".format(filts[m0],filts[m1]),fontsize=12)
        for tick in axes[indeces[0],indeces[1]].xaxis.get_major_ticks():
                tick.label.set_fontsize(8)
        # axes[indeces[0],indeces[1]].text(0.05, 0.9, \
        #     "RMS: Pull {:5.2f} Color {:6.3f}".format(numpy.std(residual/colorerror_),numpy.std(residual)), horizontalalignment='left', verticalalignment='center', \
        #     transform=axes[indeces[0],indeces[1]].transAxes,fontsize=8)
        i=i+1
fig.subplots_adjust(hspace=.4, wspace=.22)
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


correction = [fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
    + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']+ fit['zeta'][:,i][:,None]*fit['x1']\
    + fit['gamma'][:,i][:,None]*fit['k']+ fit['rho1'][:,i][:,None]*fit['R'] + (fit['ev_sig']*fit['ev'][:,i])[:,None]* fit['mag_int_raw'] \
    for i in xrange(5)]


correction = numpy.array(correction)
correction = correction - correction[2,:,:]
# correction_median = numpy.median(correction,axis=1)
from matplotlib.ticker import NullFormatter
plt.clf()
cind=[0,1,3,4]
cname = ['U','B','R','I']
mpl.rcParams['font.size'] = 12
fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'width_ratios':[.75,.25]})
for i in xrange(4):
    (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
    err = numpy.sqrt(color_cov[:,i,i] + ((ymax-ymin)/2)**2)
    #axes[i,0].errorbar(y,y-color_obs[:,i],xerr=[((ymax-ymin)/2),((ymax-ymin)/2)], yerr=[err,err],fmt='.',alpha=0.4)
    #axes[i,0].errorbar(color_obs[:,i],y-color_obs[:,i],xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[err,err],fmt='.',alpha=0.4)
    # axes[i,0].errorbar(color_obs[:,i],y,xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[((ymax-ymin)/2),((ymax-ymin)/2)],fmt='.',alpha=0.4)
    axes[i,0].errorbar(y,color_obs[:,i],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[((ymax-ymin)/2),((ymax-ymin)/2)],fmt='.',alpha=0.2)
    axes[i,0].scatter(y,color_obs[:,i],alpha=0.8,s=1)

    # axes[i,0].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]-y,xerr=[y-ymin,ymax-y], yerr=[err,err],fmt='.',alpha=0.4)
    # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.5)

 
    # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i],xerr=[y-ymin,ymax-y], yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])],fmt='.',alpha=0.4)
    # miny = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).min()
    # maxy = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).max()
    # axes[i].plot([miny,maxy],[miny,maxy])
    lims = [color_obs[:,i].min()-0.02, color_obs[:,i].max()+0.02]
    axes[i,0].plot(lims,lims,alpha=0.5)
    axes[i,0].set_xlabel(r'$\hat{{{0}}} - \hat{{V}} +(\gamma^0_\hat{{{0}}}-\gamma^0_\hat{{V}}) g_0+ (\gamma^1_\hat{{{0}}}-\gamma^1_\hat{{V}})  g_1+ (\zeta_\hat{{{0}}}-\zeta_\hat{{V}})x_1+ (\phi_\hat{{{0}}}-\phi_\hat{{V}})p$'.format(cname[i]))
    axes[i,0].set_ylabel(r'$(\hat{{{0}}}_o-\hat{{V}}_o)$'.format(cname[i]))

    axes[i,0].tick_params(axis='both', which='major', labelsize=9)
    lname = r'$(\hat{{{0}}}-\hat{{V}})$ Difference'.format(cname[i])
    # axes[i,0].set_ylabel(lname)
    # axes[i,1].hist(y-color_obs[:,i], orientation='horizontal')
    axes[i,1].hist(color_obs[:,i]-y,bins=numpy.arange(-.16,.16001,0.01))
    axes[i,1].set_xlabel(lname)
    axes[i,1].locator_params(axis='x', nbins=3)
    axes[i,1].set_xlim((-.11,.11))

    axes[i,1].xaxis.set_ticks(numpy.arange(-.10,.1001,0.10))
    axes[i,1].tick_params(axis='both', which='major', labelsize=9)
    # axes[i,1].set_ylim(axes[i,0].get_ylim())
    # axes[i,1].yaxis.set_major_formatter(NullFormatter())
    # axes[i,1].xaxis.set_major_formatter(NullFormatter())
    # axes[i].axhline(y=0,linestyle=':')
fig.subplots_adjust(hspace=.4, wspace=.18)
fig.set_size_inches(8,11)
# plt.tight_layout()
filename = 'output_fix3_x1/residual.pdf'
pp = PdfPages(filename)
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.clf()



# correction = numpy.array(correction)
# # correction = correction - correction[2,:,:]
# # correction_median = numpy.median(correction,axis=1)
# cind=[0,1,2,3,4]
# cname = ['U','B','V','R','I']
# mpl.rcParams['font.size'] = 14
# import f99_band
# A_X = f99_band.A_X(r_v=2.44,ebv=0.2/2.44)
# A_X26=A_X/(A_X[1]-A_X[2])
# A_X = f99_band.A_X(r_v=2.1,ebv=0.2/2.1)
# A_X21=A_X/(A_X[1]-A_X[2])

# fig, axes = plt.subplots(nrows=5)
# for i in xrange(5):
#     (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
#     err = numpy.sqrt(color_cov[:,1,1] + ((ymax-ymin)/2)**2)
#     #axes[i,0].errorbar(y,y-color_obs[:,i],xerr=[((ymax-ymin)/2),((ymax-ymin)/2)], yerr=[err,err],fmt='.',alpha=0.4)
#     #axes[i,0].errorbar(color_obs[:,i],y-color_obs[:,i],xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[err,err],fmt='.',alpha=0.4)
#     # axes[i,0].errorbar(color_obs[:,i],y,xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[((ymax-ymin)/2),((ymax-ymin)/2)],fmt='.',alpha=0.4)
#     axes[i].errorbar(mag_obs[:,1]-mag_obs[:,2],mag_renorm[:,i]-y,yerr=[numpy.sqrt(color_cov[:,1,1]),numpy.sqrt(color_cov[:,1,1])], xerr=[err,err],fmt='.',alpha=0.2)
#     axes[i].scatter(mag_obs[:,1]-mag_obs[:,2],mag_renorm[:,i]-y,alpha=0.8,s=1)
#     axes[i].plot(numpy.array([-0.2,0.5]),A_X26[i]*(numpy.array([-0.2,0.5])+0.025),label=r'$R_V=2.44$')
#     axes[i].plot(numpy.array([-0.2,0.5]),A_X21[i]*(numpy.array([-0.2,0.5])+0.025),label=r'$R_V=2.1$')
#     # axes[i,0].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]-y,xerr=[y-ymin,ymax-y], yerr=[err,err],fmt='.',alpha=0.4)
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.5)

 
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i],xerr=[y-ymin,ymax-y], yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])],fmt='.',alpha=0.4)

#     axes[i].set_xlabel(r'$(B_o-V_o)$')
#     axes[i].set_ylabel(r'${0}_o$'.format(cname[i]))
#     lname = r'$\Delta({0}-V)$'.format(cname[i])
#     axes[i].legend(prop={'size': 8},loc=4)

# fig.subplots_adjust(hspace=.4, wspace=.18)
# fig.set_size_inches(8,11)
# # plt.tight_layout()
# filename = 'output_fix3_x1/residual_temp.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()

# fig, axes = plt.subplots(nrows=5)
# for i in xrange(5):
#     (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
#     err = numpy.sqrt(color_cov[:,1,1] + ((ymax-ymin)/2)**2)
#     #axes[i,0].errorbar(y,y-color_obs[:,i],xerr=[((ymax-ymin)/2),((ymax-ymin)/2)], yerr=[err,err],fmt='.',alpha=0.4)
#     #axes[i,0].errorbar(color_obs[:,i],y-color_obs[:,i],xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[err,err],fmt='.',alpha=0.4)
#     # axes[i,0].errorbar(color_obs[:,i],y,xerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], yerr=[((ymax-ymin)/2),((ymax-ymin)/2)],fmt='.',alpha=0.4)
#     axes[i].errorbar(mag_obs[:,1]-mag_obs[:,2],mag_renorm[:,i],yerr=[numpy.sqrt(color_cov[:,1,1]),numpy.sqrt(color_cov[:,1,1])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.2)
#     axes[i].scatter(mag_obs[:,1]-mag_obs[:,2],mag_renorm[:,i],alpha=0.8,s=1)
#     axes[i].plot(numpy.array([-0.2,0.5]),A_X26[i]*(numpy.array([-0.2,0.5])+0.025),label=r'$R_V=2.44$')
#     axes[i].plot(numpy.array([-0.2,0.5]),A_X21[i]*(numpy.array([-0.2,0.5])+0.025),label=r'$R_V=2.1$')
#     # axes[i,0].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]-y,xerr=[y-ymin,ymax-y], yerr=[err,err],fmt='.',alpha=0.4)
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.5)

 
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i],xerr=[y-ymin,ymax-y], yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])],fmt='.',alpha=0.4)

#     axes[i].set_xlabel(r'$(B_o-V_o)$')
#     axes[i].set_ylabel(r'${0}_o$'.format(cname[i]))
#     lname = r'$\Delta({0}-V)$'.format(cname[i])
#     axes[i].legend(prop={'size': 8},loc=4)

# fig.subplots_adjust(hspace=.4, wspace=.18)
# fig.set_size_inches(8,11)
# # plt.tight_layout()
# filename = 'output_fix3_x1/residual_temp2.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()



# mpl.rcParams['font.size'] = 18
# cind=[0,1,3,4]
# cname = ['U','B','R','I']
# fig, axes = plt.subplots(nrows=4)
# for i in xrange(4):
#     (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
#     err = numpy.sqrt(color_cov[:,i,i] + ((ymax-ymin)/2)**2)
#     axes[i].errorbar(x1,color_obs[:,i]-y,xerr=[x1_err,x1_err], yerr=[err,err],fmt='.',alpha=0.4)
#     axes[i].set_xlabel(r'$X_1$'.format(cname[i]))
#     lname = r'$({0}_o-V_o) - ({0}-V)_{{model}}$'.format(cname[i])
#     axes[i].set_ylabel(lname)
#     axes[i].axhline(y=0,linestyle=':')
# fig.subplots_adjust(hspace=.3)
# fig.set_size_inches(8,11)
# filename = 'output_fix3_x1/residualx1.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


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

mpl.rcParams['font.size'] = 28
fig, axes = plt.subplots()
bins = numpy.arange(-0.2,0.8001,0.02)
plt.hist(fit['Delta'].flatten(),bins,label='ideogram',normed=True,alpha=0.5)
plt.hist(numpy.median(fit['Delta'],axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
plt.legend()
plt.xlabel(r'$\Delta$')
pp = PdfPages('output_fix3_x1/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

fig, axes = plt.subplots()
# import matplotlib.pyplot
# matplotlib.pyplot.rcdefaults()
# dmean = fit['Delta'].mean()
x, xmin, xmax = numpy.percentile(fit['Delta'],(50, 50-34,50+34),axis=0)
dum = x>0.39
print zhelio[dum],x[dum],snname[dum]

# fig, axes = plt.subplots()
# plt.scatter(x,xmin)
# plt.xlabel(r'$\Delta$')
# plt.ylabel(r'$-1\sigma \Delta$')
# pp = PdfPages('output_fix3_x1/Delta1sigvsDelta.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# fig, axes = plt.subplots()
# plt.scatter(x,(x-xmin)/(xmax-x))
# plt.xlabel(r'$\Delta$')
# plt.ylabel(r'$-1\sigma / +1\sigma$')
# pp = PdfPages('output_fix3_x1/sigmadiffvsDelta.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# x, xmin, xmax = numpy.percentile(fit['Delta']-fit['Delta'][:,0][:,None],(50, 50-34,50+34),axis=0)

# fig, axes = plt.subplots()
# plt.scatter(x,(x-xmin)/(xmax-x))
# plt.xlabel(r'$\delta \Delta$')
# plt.ylabel(r'$-1\sigma / +1\sigma$')
# pp = PdfPages('output_fix3_x1/deltasigmadiffvsDelta.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# fig, axes = plt.subplots()
# plt.errorbar(zhelio,x,yerr=(x-xmin,xmax-x),fmt='o')
# plt.ylabel(r'$\delta \Delta$')
# plt.xlabel(r'$z_{\odot}$')
# pp = PdfPages('output_fix3_x1/deltaDelta_vs_z.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# fig, axes = plt.subplots()
# plt.hist((fit['Delta']-fit['Delta'][:,0][:,None]).flatten(),bins=20)
# plt.ylabel(r'$\delta \Delta$')
# pp = PdfPages('output_fix3_x1/deltaDelta_hist.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()
print numpy.std(numpy.median((fit['Delta']-fit['Delta'][:,0][:,None]),axis=0))

fig, axes = plt.subplots()
bins = numpy.arange(-0.2,0.8001,0.02)
plt.hist((fit['Delta']-fit['Delta'][:,0][:,None]).flatten(),bins,label='posterior stack',normed=True,alpha=0.5)
plt.hist(numpy.median((fit['Delta']-fit['Delta'][:,0][:,None]),axis=0),bins,label='median',normed=True,alpha=0.5,width=0.01)
plt.legend(fontsize=20)
plt.xlabel(r'$\Delta-\Delta|_0$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
pp = PdfPages('output_fix3_x1/deltaDelta_hist.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()



x, xmin, xmax = numpy.percentile(fit['Delta']-fit['Delta'][:,0][:,None],(50,50-34,50+34),axis=0)
plt.errorbar(zhelio,x,yerr=(x-xmin,xmax-x),fmt='o')
plt.ylabel(r'$\Delta$')
plt.xlabel(r'$z_{\odot}$')
pp = PdfPages('output_fix3_x1/Delta_vs_z.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.errorbar(color_obs[:,0],x,yerr=(x-xmin,xmax-x),fmt='o')
plt.ylabel(r'$\Delta$')
plt.xlabel(r'$U-V$')
pp = PdfPages('output_fix3_x1/Delta_vs_u-v.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.errorbar(color_obs[:,3],x,yerr=(x-xmin,xmax-x),fmt='o')
plt.ylabel(r'$\Delta$')
plt.xlabel(r'$I-V$')
pp = PdfPages('output_fix3_x1/Delta_vs_i-v.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


fig, axes = plt.subplots()
plt.hist(numpy.median(fit['Delta'],axis=0),bins=20)
plt.xlabel(r'$\Delta$')
plt.ylabel(r'Number per bin')
pp = PdfPages('output_fix3_x1/Delta_med_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['k'].flatten(),normed=True,bins=20)
plt.title(r'$k$')
pp = PdfPages('output_fix3_x1/k_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['R'].flatten(),normed=True,bins=20)
plt.title(r'$D$')
pp = PdfPages('output_fix3_x1/D_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1]-mag_obs[:,2])
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$B_o-V_o$')
pp = PdfPages('output_fix3_x1/ebvvsobs.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1],color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29. + 3.96*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$B_o$')
pp = PdfPages('output_fix3_x1/g1.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,0],color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$U_o$')
pp = PdfPages('output_fix3_x1/g1u.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0), \
    mag_obs[:,0]- numpy.median(fit['alpha'][:,0][:,None]*fit['EW'][:,:,0]) - numpy.median(fit['beta'][:,0][:,None]*fit['EW'][:,:,1])\
    ,color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$U_o - \alpha_0 EW_{Ca} - \beta_0 EW_{Si}$')
pp = PdfPages('output_fix3_x1/g2uc.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.scatter(numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],axis=0),mag_obs[:,1]-3.96*numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),color='blue')
plt.plot(numpy.array([-0.015,0.08]), -29 - 3.46*numpy.array([-0.015,0.08]))
plt.xlabel(r'$E_\delta(B-V)$')
plt.ylabel(r'$B_o - \gamma_1 E_\gamma(B-V)$')
pp = PdfPages('output_fix3_x1/g2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



# plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']],normed=True,bins=20,
#     label=[r'$E_\gamma(B-V)$',r'$E_\delta(B-V)$'],range=(-0.1,0.4))
# plt.xlabel(r'$E(B-V)$',fontsize=20)
# plt.legend()
# pp = PdfPages('output_fix3_x1/ebv.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

dustebv = (fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']
extebv = (fit['ev_sig']*(fit['ev'][:,1]-fit['ev'][:,2]))[:,None]*fit['mag_int_raw']
dustebv = dustebv-dustebv[:,0][:,None]
extebv = extebv-extebv[:,0][:,None]


plt.clf()
plt.figure(figsize=(8.0, 6.0))
plt.hist([dustebv, extebv],normed=True,bins=50,
  label=[r'$E_{\gamma}(\hat{B}-\hat{V})$',r'$E_p(\hat{B}-\hat{V})$'])
# plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] + 
#   (fit['gamma1'][:,1]-fit['gamma1'][:,2])[:,None]*fit['k1'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']],normed=True,bins=20,
#   label=[r'$E(B-V)$',r'$E_\delta(B-V)$'])
plt.xlabel(r'$E(\hat{B}-\hat{V})-E(\hat{B}-\hat{V})|_0$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=20)

pp = PdfPages('output_fix3_x1/ebv.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


plt.hist( extebv.flatten(),normed=True,bins=bins,alpha=0.5,
  label='ideogram')
plt.hist( numpy.median(extebv,axis=0),normed=True,bins=bins,alpha=0.5,width=0.005,
  label='median')
plt.xlabel(r'$E_\delta(\hat{B}-\hat{V})-E_\delta(\hat{B}-\hat{V})|_0$',fontsize=20)
plt.xlim((-.1,.12))
plt.legend()
pp = PdfPages('output_fix3_x1/ebv_phi.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


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
# pp = PdfPages("output_fix3_x1/Rveff.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (y, ymin, ymax) = numpy.percentile((au-av)/(ab-av),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.ylim((1.8,2.3))
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$R_{eff,U} - R_{eff,V} = E_{eff}(U-V) /E_{eff}(U-V)$')
# pp = PdfPages("output_fix3_x1/Rveff2.pdf")
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
# (y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$R_{eff, V} = A_{eff, V}/E_{eff}(B-V)$')
# pp = PdfPages('output_fix3_x1/Rveff2.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# (x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
# (y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
# plt.errorbar(x,y/(ymax-ymin)*2,xerr=[x-xmin,xmax-x],fmt='o')
# plt.xlim((-0.04,0.5))
# plt.xlabel(r'$E_{eff}(B-V)$')
# plt.ylabel(r'$S/N(R_{eff, V})$')
# pp = PdfPages('output_fix3_x1/Rveff3.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.errorbar(x,y-3.96*x,xerr=[x-xmin,xmax-x],yerr=[(y-ymin),(ymax-y)],fmt='o')
# plt.xlabel(r'$E_{eff}(B-V)+const$')
# plt.ylabel(r'$\Delta A_B  = A_{eff, B}-3.96 (E_{eff}(B-V)+const)$')
# plt.ylim((-0.8,0.1))
# pp = PdfPages('output_fix3_x1/Rveffres.pdf')
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
# pp = PdfPages('output_fix3_x1/ebvebv.pdf')
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
# pp = PdfPages('output_fix3_x1/Rveff.pdf')
# plt.xlabel(r'$R_{V,eff}$')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median(ab-av,axis=0), numpy.median(rv,axis=0))
# plt.xlabel(r'$E^o(B-V)$')
# plt.ylabel(r'$R_{V,eff}$')
# pp = PdfPages('output_fix3_x1/EVRveff.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
pp = PdfPages('output_fix3_x1/c_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages('output_fix3_x1/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['alpha']-fit['alpha'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$"])
# pp = PdfPages('output_fix3_x1/alpham4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()



figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages('output_fix3_x1/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['eta'],labels=[r"${\eta}_0$",r"${\eta}_1$",r"${\eta}_2$",r"${\eta}_3$",r"${\eta}_4$"])
pp = PdfPages('output_fix3_x1/eta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['beta']-fit['beta'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$"])
# pp = PdfPages('output_fix3_x1/betam4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


figure = corner.corner(fit['gamma'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_2$",r"${\gamma}_3$",r"${\gamma}_4$"])
pp = PdfPages('output_fix3_x1/gamma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['rho1'],labels=[r"${\rho}_{10}$",r"${\rho}_{11}$",r"${\rho}_{12}$",r"${\rho}_{13}$",r"${\rho}_{14}$"])
pp = PdfPages('output_fix3_x1/rho_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


with PdfPages('output_fix3_x1/multipage_pdf.pdf') as pdf:

    lineobjects = plt.plot(fit['lp__'][::10])
    plt.title(r'log p')
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
    pp = PdfPages('output_fix3_x1/rx_corner.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()

    mega = fit['rho1']/((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None]

    figure = corner.corner(mega,labels=[r"$R_{\delta U}$",r"$R_{\delta B}$",r"$R_{\delta V}$",r"$R_{\delta R}$",r"$R_{\delta I}$"], \
        range=[[-8.,0.5] for x in xrange(5)])
    pp = PdfPages('output_fix3_x1/rxdelta_corner.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()


#from manu
efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])
# [3701, 4601, 5744, 6948, 8403]


filts = [r'$\hat{U}$',r'$\hat{B}$',r'$\hat{V}$',r'$\hat{R}$',r'$\hat{I}$']


labels = [r'$\hat{U}$',r'$\hat{B}$',r'$\hat{V}$',r'$\hat{R}$',r'$\hat{I}$']
from matplotlib.ticker import FuncFormatter, MaxNLocator
def format_fn2(tick_val, tick_pos):
    if int(tick_val) in numpy.arange(5):
        return labels[int(tick_val)]
    else:
        return ''

(y, ymin, ymax) = numpy.percentile(fit['ev']/fit['ev'][:,2][:,None],(50,50-34,50+34),axis=0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(numpy.arange(5),y-1,yerr=[y-ymin,ymax-y],fmt='o')
ax.xaxis.set_major_formatter(FuncFormatter(format_fn2))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.axhline(0,linestyle=':')
ax.set_xlabel(r'Band $X$')
ax.set_xlim((-0.5,4.5))
ax.set_ylabel(r'$\frac{\phi_X}{\phi_{\hat{V}}}-1$')
ax.set_ylim((-1.6,0.1))
pp = PdfPages('output_fix3_x1/phiratio.pdf')
plt.savefig(pp,format='pdf',bbox_inches='tight')
pp.close()
plt.close()


mega = numpy.array([EW_obs[:,0],EW_obs[:,0],sivel, x1,c ])


mega = numpy.transpose(mega)
c=ChainConsumer()
c.add_chain(mega, \
    parameters= [r"$EW_{Ca}$",r"$EW_{Si}$",r"$\lambda_{Si}$",r"$x_1$",r"$C$"],name='This Article')
d= [x.split() for x in open("/Users/akim/project/Pantheon/data_fitres/Ancillary_C11.FITRES").readlines()]
panth = []
for dl in d[69:]:
    panth.append([float(dl[7]), float(dl[20]), float(dl[22])])

panth=numpy.array(panth)
goodz  = numpy.logical_and(panth[:,0] > 0.03, panth[:,0] < 0.08)
c.add_chain(panth[goodz,1:],parameters= [r"$x_1$",r"$C$"],name=r"Pantheon $0.03<z<0.08$")
c.configure(bins=10)

fig = c.plotter.plot(figsize="column", truth=numpy.zeros(5))
# fig = c.plotter.plot_distributions()

for ax in fig.axes:
    ax.xaxis.set_tick_params(labelsize=6)
    ax.xaxis.label.set_size(7)
    ax.yaxis.set_tick_params(labelsize=6)
    ax.yaxis.label.set_size(7)
fig.savefig('output_fix3_x1/perobject_input.pdf',bbox_inches='tight')


mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), fit['x1'].flatten(),\
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']).flatten(), \
    ((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']).flatten(),((fit['ev_sig']*fit['ev'][:,2])[:,None]* fit['mag_int_raw']).flatten()])


mega = numpy.transpose(mega)
mega=mega[::50,:]
mega[:,1] = mega[:,1]+EW_mn[0]
mega[:,2] = mega[:,2]+EW_mn[1]
mega[:,3] = mega[:,3]+sivel_mn



c=ChainConsumer()
c.add_chain(mega, \
    parameters= [r"$\Delta$",r"$EW_{Ca}$",r"$EW_{Si}$",r"$\lambda_{Si}$",r"$x_1$",r"$E_{\gamma^0}(B-V)$",r"$E_{\gamma^1}(B-V)$",r"$A_{p,V}$"])

fig = c.plotter.plot(figsize="column", truth=numpy.zeros(5))

for ax in fig.axes:
    ax.xaxis.set_tick_params(labelsize=4)
    ax.xaxis.label.set_size(5)
    ax.yaxis.set_tick_params(labelsize=4)
    ax.yaxis.label.set_size(5)
fig.savefig('output_fix3_x1/perobject_corner.pdf',bbox_inches='tight')

cornered

# figure = corner.corner(mega,labels=[r"$\Delta$",r"$EW_{Ca}$",r"$EW_{Si}$",r"$\lambda_{Si}$",r"$x_1$",r"$E_{\gamma^0}(B-V)$",r"$E_{\gamma^1}(B-V)$",r"$A_{p,V}$"],range=numpy.zeros(8)+.9995,label_kwargs={'fontsize':22})
# for ax in figure.get_axes():
#     for tick in ax.xaxis.get_major_ticks():
#         tick.label.set_fontsize(14) 
#     for tick in ax.yaxis.get_major_ticks():
#         tick.label.set_fontsize(14) 
# pp = PdfPages('output_fix3_x1/perobject_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
plt.close()
plt.clf()


nlinks = fit['gamma'].shape[0]
mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['zeta'],fit['gamma'],fit['rho1'],fit['ev_sig'][:,None]*fit['ev']])
mega = numpy.transpose(mega)

cname=[r'\hat{U}',r'\hat{B}',r'\hat{V}',r'\hat{R}',r'\hat{I}']
for index in xrange(5):
    c=ChainConsumer()
    c.add_chain(mega[index,:,:],parameters=[r"$c_{{{}}}$".format(cname[index]), r"$\alpha_{{{}}}$".format(cname[index]),\
                    r"$\beta_{{{}}}$".format(cname[index]),r"$\eta_{{{}}}$".format(cname[index]),r"$\zeta_{{{}}}$".format(cname[index]),r"$\gamma^0_{{{}}}$".format(cname[index]),\
                    r"$\gamma^1_{{{}}}$".format(cname[index]),r"$\sigma_p \phi_{{{}}}$".format(cname[index])])
    fig = c.plotter.plot( figsize="column", truth=[None,0,0,0,0,0,0,0])

    for ax in fig.axes:
        ax.xaxis.set_tick_params(labelsize=7)
        ax.xaxis.label.set_size(9)
        ax.yaxis.set_tick_params(labelsize=7)
        ax.yaxis.label.set_size(9)
    fig.savefig(filename='output_fix3_x1/coeff{}.pdf'.format(index),bbox_inches='tight')

c=ChainConsumer()
c.add_chain(fit['ev_sig'][:,None]*fit['ev']*numpy.sign(fit['ev'][:,2])[:,None],parameters=[r"$\sigma_p \phi_{\hat{U}}$",r"$\sigma_p \phi_{\hat{B}}$",r"$\sigma_p \phi_{\hat{V}}$",r"$\sigma_p \phi_{\hat{R}}$",r"$\sigma_p \phi_{\hat{I}}$"])
c.plotter.plot(filename='output_fix3_x1/sigev.pdf',figsize="column", truth=numpy.zeros(5))




