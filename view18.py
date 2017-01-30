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

ext=''

f = open('temp18AV'+ext+'.pkl','rb')
(fit,_) = pickle.load(f)

# for key in fit.keys():
#     print key, fit[key].min(), fit[key].max()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

sivel,sivel_err,x1,x1_err, _, _ = sivel.sivel(data)
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
EW_obs=EW_obs[use]
mag_obs=mag_obs[use]
EW_cov= EW_cov[use]
mag_cov=mag_cov[use]

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

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))


figure = corner.corner(numpy.concatenate((EW_obs,sivel[:,None]),axis=1),labels=[r"$EW_{Ca,o}$",r"$EW_{Si,o}$",r"$\lambda_{Si,o}$"])
for ax in figure.get_axes():
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(18) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        
pp = PdfPages('output18'+ext+'/feature_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# colors = fit['c']+mag_mn[None,:]
# colors  = colors - colors[:,2][:,None]
# colors = numpy.delete(colors, 2,axis=1)

# figure = corner.corner(colors,labels=[r"$U_0-V_0$",r"$B_0-V_0$",r"$R_0-V_0$",r"$I_0-V_0$"])
# pp = PdfPages('output18'+ext+'/col_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


(y, ymin, ymax) = numpy.percentile(fit['EW'][:,:,1],(50,50-34,50+34),axis=0)
plt.errorbar(x1, y, xerr=[x1_err,x1_err],yerr=[y-ymin,ymax-ymin],fmt='o')
plt.xlabel(r'$X_1$')
plt.ylabel(r'$EW_{Si}$')
pp = PdfPages('output18'+ext+'/x1si.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()





# correction = [fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
#     + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] + fit['eta'][:,i][:,None]*fit['sivel']\
#     + fit['gamma'][:,i][:,None]*fit['k']+ fit['rho1'][:,i][:,None]*fit['R'] \
#     for i in xrange(5)]

# correction = numpy.array(correction)
# correction = correction - correction[2,:,:]
# # correction_median = numpy.median(correction,axis=1)

# cind=[0,1,3,4]
# cname = ['U','B','R','I']
# fig, axes = plt.subplots(nrows=4)
# for i in xrange(4):
#     (y, ymin, ymax) = numpy.percentile(correction[cind[i],:,:],(50,50-34,50+34),axis=0)
#     # err = numpy.sqrt(color_cov[:,i,i] + ((ymax-ymin)/2)**2)
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]-y,xerr=[y-ymin,ymax-y], yerr=[err,err],fmt='.',alpha=0.4)
#     axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2],yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])], xerr=[(ymax-ymin)/2,(ymax-ymin)/2],fmt='.',alpha=0.5)
#     # axes[i].errorbar(y+mag_mn[cind[i]]-mag_mn[2],color_obs[:,i],xerr=[y-ymin,ymax-y], yerr=[numpy.sqrt(color_cov[:,i,i]),numpy.sqrt(color_cov[:,i,i])],fmt='.',alpha=0.4)
#     miny = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).min()
#     maxy = (color_obs[:,i]+mag_mn[cind[i]]-mag_mn[2]).max()
#     axes[i].plot([miny,maxy],[miny,maxy])
#     axes[i].set_ylabel(r'$({0}_o-V_o)$'.format(cname[i]))
#     lname = r'$({0}-V)_{{model}}$'.format(cname[i])
#     axes[i].set_xlabel(lname)
#     # axes[i].axhline(y=0,linestyle=':')
# fig.subplots_adjust(hspace=.3)
# fig.set_size_inches(8,11)
# filename = 'output18'+ext+'/residual.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()


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
# filename = 'output18'+ext+'/residualx1.pdf'
# pp = PdfPages(filename)
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()



# print 'Delta', numpy.median(fit['Delta'][:,outlier[0]],axis=0)
# print 'k', numpy.median(fit['k'][:,outlier[0]],axis=0)
# print 'EW0', numpy.median(fit['EW'][:,outlier[0],0],axis=0), EW_renorm[outlier[0],0] 
# print 'EW1', numpy.median(fit['EW'][:,outlier[0],1],axis=0), EW_renorm[outlier[0],1]
fig, axes = plt.subplots()
plt.hist(fit['Delta'].flatten(),normed=True,bins=20,color='yellow')
plt.xlabel(r'$\Delta$')
pp = PdfPages('output18'+ext+'/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['AV'].flatten(),normed=True,bins=20)
plt.xlabel(r'$A_V$')
pp = PdfPages('output18'+ext+'/AV_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(numpy.median(fit['AV'],axis=0),normed=True,bins=20)
plt.xlabel(r'$A_V$')
pp = PdfPages('output18'+ext+'/AV_mode_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['RVinv'].flatten(),normed=True,bins=20)
plt.xlabel(r'$R_V^{-1}$')
pp = PdfPages('output18'+ext+'/RVinv_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(1./fit['RVinv'].flatten(),normed=True,bins=20)
plt.xlabel(r'$R_V$')
pp = PdfPages('output18'+ext+'/RV_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# plt.hist((fit['AV']*fit['RVinv']).flatten(),normed=True,bins=20)
# plt.xlabel(r'$E(B-V)$')
# plt.xlim((-1,5))
# pp = PdfPages('output18'+ext+'/EBV_hist.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
# pp = PdfPages('output18'+ext+'/c_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages('output18'+ext+'/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages('output18'+ext+'/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['eta'],labels=[r"${\eta}_0$",r"${\eta}_1$",r"${\eta}_2$",r"${\eta}_3$",r"${\eta}_4$"])
pp = PdfPages('output18'+ext+'/eta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# figure = corner.corner(fit['deltaAV'],labels=[r"${\delta}_{AV,0}$",r"${\delta}_{AV,1}$",r"${\delta}_{AV,2}$",r"${\delta}_{AV,3}$",r"${\delta}_{AV,4}$"])
# pp = PdfPages('output18'+ext+'/deltaAV_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


# mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages('output18'+ext+'/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


with PdfPages('output18'+ext+'/multipage_pdf.pdf') as pdf:

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



    # plt.plot(fit['Delta_scale'][::10])
    # plt.title(r'$\Delta$ scale')
    # pdf.savefig()
    # plt.close()
    
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
 

#from manu
efflam = numpy.array([ 3693.16777627,  4369.37505509,  5287.48667023,  6319.19906153,7610.89305298])
# [3701, 4601, 5744, 6948, 8403]
rc('text', usetex=True)





# for index in xrange(5):
#     figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index), r"$\alpha_{}$".format(index),\
#                     r"$\beta_{}$".format(index),r"$\eta_{}$".format(index), r"$\sigma_{}$".format(index)],label_kwargs={'fontsize':20})
#     figure.suptitle(filts[index],fontsize=28)

#     for ax in figure.get_axes():
#         for tick in ax.xaxis.get_major_ticks():
#             tick.label.set_fontsize(18'+ext+') 
#         for tick in ax.yaxis.get_major_ticks():
#             tick.label.set_fontsize(18) 

#     pp = PdfPages('output18/coeff{}.pdf'.format(index))
#     plt.savefig(pp,format='pdf')
#     pp.close()
#     plt.close()



mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), \
    fit['AV'].flatten()])

mega = numpy.transpose(mega)
mega=mega[::50,:]

figure = corner.corner(mega,labels=[r"$\Delta$",r"$EW_{Ca}$",r"$EW_{Si}$",r"$\lambda_{Si}$",r"$A_V$"],range=numpy.zeros(5)+1.)
for ax in figure.get_axes():
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(18) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(18) 
pp = PdfPages('output18'+ext+'/perobject_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


