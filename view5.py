#!/usr/bin/env python
import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo

f = open('temp5.pkl','rb')
(fit,_) = pickle.load(f)

for key in fit.keys():
    print key, fit[key].min(), fit[key].max()

# flip = fit['rho0'][:,4] <0
# for i in xrange(5):
#     fit['rho0'][flip,i] *= -1

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

# dum = numpy.zeros((5,5))
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

# dum/= len(fit['L_Omega'])
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

# trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
# trans = numpy.array(trans)
# color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in color_cov])

# correction = [  fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
#     + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] \
#     + fit['gamma'][:,i][:,None]*fit['k'] \
#     + fit['Delta'] \
#     + fit['rho0'][:,i][:,None]*fit['R'] \
#     + fit['rho1'][:,i][:,None]*fit['R']**2 \
#     - fit['rho1'][:,i][:,None] \
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

print fit['Delta'].flatten().std()
plt.hist(fit['Delta'].flatten(),normed=True,bins=20)
plt.title(r'$\Delta$')
pp = PdfPages('output5/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['R'].flatten(),normed=True,bins=20)
plt.title(r'$\delta$')
pp = PdfPages('output5/D_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['R'].flatten(),normed=True,bins=20)
plt.title(r'$k$')
pp = PdfPages('output5/k_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


# plt.hist(fit['R_scale'].flatten(),normed=True,bins=20)
# plt.title(r'R_scale')
# pp = PdfPages('output5/Rscale_hist.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
pp = PdfPages('output5/c_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages('output5/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages('output5/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['gamma'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_2$",r"${\gamma}_3$",r"${\gamma}_4$"])
pp = PdfPages('output5/gamma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages('output5/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# figure = corner.corner(fit['rho0'],labels=[r"${\delta}_{00}$",r"${\delta}_{01}$",r"${\delta}_{02}$",r"${\delta}_{03}$",r"${\delta}_{04}$"])
# pp = PdfPages('output5/delta0_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['rho1'],labels=[r"${\delta}_{10}$",r"${\delta}_{11}$",r"${\delta}_{12}$",r"${\delta}_{13}$",r"${\delta}_{14}$"])
pp = PdfPages('output5/delta1_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
# pp = PdfPages('output5/c_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

mega = fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None]
figure = corner.corner(mega,labels=[r"$R_0$",r"$R_1$",r"$R_2$",r"$R_3$",r"$R_4$"])
pp = PdfPages('output5/rx_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

mega = fit['rho1']/((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None]
figure = corner.corner(mega,labels=[r"$R_{\rho 0}$",r"$R_{\rho 1}$",r"$R_{\rho 2}$",r"$R_{\rho 3}$",r"$R_{\rho 4}$"])
pp = PdfPages('output5/rx2_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

with PdfPages('output5/multipage_pdf.pdf') as pdf:

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
    # plt.title(r'$\Delta$_scale')
    # pdf.savefig()
    # plt.close()
    # plt.plot(fit['Delta_scale'][::10])
    # plt.title(r'$\Delta scale$')
    # pdf.savefig()
    # plt.close()

    # lineobjects = plt.plot(fit['c'][::10,:],label=['U','B','V','R','I'])
    # plt.title('c')
    # plt.legend(iter(lineobjects),('U','B','V','R','I'))
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
    
    lineobjects = plt.plot(fit['gamma'][::10],label=['U','B','V','R','I'])
    plt.title(r'$\gamma$')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()
 
    # lineobjects = plt.plot(fit['rho0'][::10],label=['U','B','V','R','I'])
    # plt.title(r'$\delta0$')
    # plt.legend(iter(lineobjects),('U','B','V','R','I'))
    # pdf.savefig()
    # plt.close()

    lineobjects = plt.plot(fit['rho1'][::10],label=['U','B','V','R','I'])
    plt.title(r'$\delta1$')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(fit['L_sigma'][::10],label=['U','B','V','R','I'])
    plt.title(r'$\sigma$')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
    pdf.savefig()
    plt.close()

    # plt.plot(fit['R'][::10,:10])
    # plt.title('R')
    # pdf.savefig()
    # plt.close()
    


nlinks = fit['gamma'].shape[0]
dum = numpy.concatenate((fit['gamma'][:,0:2], numpy.random.normal(1.,0.01,size=(nlinks,1)),fit['gamma'][:,3:]),axis=1)
mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],dum,fit['rho1'],fit['L_sigma']])
mega = numpy.transpose(mega)

for index in xrange(5):
    figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index), r"$\alpha_{}$".format(index),\
                    r"$\beta_{}$".format(index),r"$\gamma_{}$".format(index),\
                    r"$\delta_{{1{}}}$".format(index),r"$\sigma_{}$".format(index)])
    pp = PdfPages('output5/coeff{}.pdf'.format(index))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()


    # filter wavelengths
r = numpy.log(8635./3276.)
edges= 3276.*numpy.exp(numpy.arange(6)*r/5.)
# edges = numpy.array([[3300., 4102], [4102, 5100], [5200, 6289], [6289, 7607], [7607, 9200]])
efflam = []
for i in xrange(5):
    efflam.append((edges[i+1]-edges[i])/2+edges[i])
    
# [3701, 4601, 5744, 6948, 8403]
rc('text', usetex=True)


rvs = [1.7,3.1,4.1]
filts = ['U','B','V','R','I']
fig, axes = plt.subplots(nrows=len(filts),sharex=True)
xerr = numpy.sqrt(EW_cov[:,0,0]) 
for i in xrange(len(filts)):
    r = numpy.array( [numpy.argmin(EW_obs[:,0]), numpy.argmax(EW_obs[:,0])])
    yerr= numpy.sqrt(mag_cov[:,i,i])
    axes[i].errorbar(EW_obs[:,0],mag_obs[:, i], \
        xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
    axes[i].set_ylabel(r'{}'.format(filts[i]))
    offset  = numpy.mean(mag_obs[:, i] - numpy.median(fit['alpha'][:,i][:,None]*fit['EW'][:,:,0],axis=0))
    axes[i].plot(EW_obs[r,0],offset+numpy.median(fit['alpha'][:,i],axis=0)*EW_renorm[r,0] \
        ,color='red',linewidth=2)
axes[len(filts)-1].set_xlabel(r'EW(Ca)')
fig.subplots_adjust(hspace=0.001)
pp = PdfPages('output5/speccamag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

fig, axes = plt.subplots(nrows=len(filts),sharex=True)
xerr = numpy.sqrt(EW_cov[:,1,1])
for i in xrange(len(filts)):
    r = numpy.array( [numpy.argmin(EW_obs[:,1]), numpy.argmax(EW_obs[:,1])])
    yerr= numpy.sqrt(mag_cov[:,i,i])
    axes[i].errorbar(EW_obs[:,1],mag_obs[:, i], \
        xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
    axes[i].set_ylabel(r'{}'.format(filts[i]))
    offset  = numpy.mean(mag_obs[:, i] - numpy.median(fit['beta'][:,i][:,None]*fit['EW'][:,:,1],axis=0))
    axes[i].plot(EW_obs[r,1],offset+numpy.median(fit['beta'][:,i],axis=0)*EW_renorm[r,1] \
        ,color='red',linewidth=2)
axes[len(filts)-1].set_xlabel(r'EW(Si)')
fig.subplots_adjust(hspace=0.001)
pp = PdfPages('output5/specsimag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()
fig, axes = plt.subplots(nrows=len(filts),sharex=True)
xerr = numpy.sqrt(mag_cov[:,1,1]+ mag_cov[:,2,2] - 2*mag_cov[:,1,2])
for i in xrange(len(filts)):
    r = numpy.array( [numpy.argmin(mag_obs[:,1]-mag_obs[:,2]), numpy.argmax(mag_obs[:,1]-mag_obs[:,2])])
    yerr= numpy.sqrt(mag_cov[:,i,i])
    axes[i].errorbar(mag_obs[:,1]-mag_obs[:,2],mag_obs[:, i], \
        xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
    axes[i].set_ylabel(r'{}'.format(filts[i]))
    slope = numpy.median(fit['gamma'][:,i] / (fit['gamma'][:,1]-fit['gamma'][:,2]),axis=0)
    offset  = numpy.mean(mag_obs[:, i] - numpy.median(slope*(mag_obs[:,1]-mag_obs[:,2])))
    axes[i].plot(mag_obs[r,1]-mag_obs[r,2],offset+slope*(mag_obs[r,1]-mag_obs[r,2]) \
        ,color='red',linewidth=2)
axes[len(filts)-1].set_xlabel(r'B-V')
fig.subplots_adjust(hspace=0.001)
pp = PdfPages('output5/colormag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


lambdas = numpy.arange(3500.,9000,100)
for rv in rvs:
    A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
    norm  = sncosmo._extinction.ccm89(numpy.array([efflam[1]]), 1., rv)-sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
    A_ = A_/norm[0]
    plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))
(y, ymin, ymax) = numpy.percentile(fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]),(50,50-34,50+34),axis=0)

plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
plt.legend()
plt.xlabel(r'Wavelength (\AA)')
plt.ylabel(r'$R_X=\frac{\gamma_X}{\gamma_1-\gamma_2}$')
pp = PdfPages('output5/ccm.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

for rv in rvs:
    A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
    norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
    A_ = A_/norm[0]
    plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))
(y, ymin, ymax) = numpy.percentile(fit['gamma']/fit['gamma'][:,2][:,None],(50,50-34,50+34),axis=0)

plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
plt.legend()
plt.xlabel(r'Wavelength (\AA)')
plt.ylabel(r'$\frac{\gamma_X}{\gamma_2}$')
pp = PdfPages('output5/ccm2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# for index in xrange(5):
#     dum = numpy.swapaxes(numpy.array([fit['alpha'][:,index]*color_std,fit['beta'][:,index]*color_std]),0,1)
#     figure = corner.corner(dum,labels=[r"$\alpha_{}$".format(index),r"$\beta_{}$".format(index)],truths=[0,0])
#     pdf.savefig()
#     plt.close()

# plt.plot(fit['ebeta_inv'])
# plt.title('ebeta_inv')
# pdf.savefig()
# plt.close()