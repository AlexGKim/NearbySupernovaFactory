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


f = open('temp11.pkl','rb')
(fit,_) = pickle.load(f)

# for key in fit.keys():
#     print key, fit[key].min(), fit[key].max()

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

EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

sivel_mn = sivel.mean()
sivel_renorm = sivel-sivel_mn


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

plt.hist(fit['Delta'].flatten(),normed=True,bins=20)
plt.xlabel(r'$\Delta$')
pp = PdfPages('output11/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['k'].flatten(),normed=True,bins=20)
plt.title(r'$k$')
pp = PdfPages('output11/k_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['R'].flatten(),normed=True,bins=20)
plt.title(r'$D$')
pp = PdfPages('output11/D_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1]-mag_obs[:,2])
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$B_o-V_o$')
pp = PdfPages('output11/ebvvsobs.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,1],color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29. + 3.96*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$B_o$')
pp = PdfPages('output11/g1.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),mag_obs[:,0],color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$U_o$')
pp = PdfPages('output11/g1u.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.scatter(numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0), \
    mag_obs[:,0]- numpy.median(fit['alpha'][:,0][:,None]*fit['EW'][:,:,0]) - numpy.median(fit['beta'][:,0][:,None]*fit['EW'][:,:,1])\
    ,color='red',label=r'$E_\gamma(B-V)$')
plt.plot(numpy.array([-0.07,0.36]), -29.2 + 4.87*numpy.array([-0.08,0.36]))
plt.xlabel(r'$E_\gamma(B-V)$')
plt.ylabel(r'$U_o - \alpha_0 EW_{Ca} - \beta_0 EW_{Si}$')
pp = PdfPages('output11/g2uc.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


plt.scatter(numpy.median((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'],axis=0),mag_obs[:,1]-3.96*numpy.median((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'],axis=0),color='blue')
plt.plot(numpy.array([-0.015,0.08]), -29 - 3.46*numpy.array([-0.015,0.08]))
plt.xlabel(r'$E_\delta(B-V)$')
plt.ylabel(r'$B_o - \gamma_1 E_\gamma(B-V)$')
pp = PdfPages('output11/g2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



plt.hist([(fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'], (fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']],normed=True,bins=20,
    label=[r'$E_\gamma(B-V)$',r'$E_\delta(B-V)$'],range=(-0.1,0.4))
plt.xlabel(r'$E(B-V)$')
plt.legend()
pp = PdfPages('output11/ebv.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

ab = fit['gamma'][:,1][:,None]*fit['k'] + fit['rho1'][:,1][:,None]*fit['R']
av = fit['gamma'][:,2][:,None]*fit['k'] + fit['rho1'][:,2][:,None]*fit['R']
rv = av/(ab-av+0.2)
(x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
(y, ymin, ymax) = numpy.percentile(av,(50,50-34,50+34),axis=0)+0.15
plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
x_=numpy.array([-0.03,0.45])
plt.plot(x_,2.96*x_)
plt.xlim((-0.04,0.5))
plt.xlabel(r'$E_{eff}(B-V)$')
plt.ylabel(r'$A_{eff, V}$')
pp = PdfPages('output11/Rveff.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
(y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
plt.xlim((-0.04,0.5))
plt.xlabel(r'$E_{eff}(B-V)$')
plt.ylabel(r'$R_{eff, V} = A_{eff, V}/E_{eff}(B-V)$')
pp = PdfPages('output11/Rveff2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(x, xmin, xmax) = numpy.percentile(ab-av,(50,50-34,50+34),axis=0) + 0.08
(y, ymin, ymax) = numpy.percentile((av+0.15)/(ab-av+0.08),(50,50-34,50+34),axis=0)
plt.errorbar(x,y/(ymax-ymin)*2,xerr=[x-xmin,xmax-x],fmt='o')
plt.xlim((-0.04,0.5))
plt.xlabel(r'$E_{eff}(B-V)$')
plt.ylabel(r'$S/N(R_{eff, V})$')
pp = PdfPages('output11/Rveff3.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.errorbar(x,y-3.96*x,xerr=[x-xmin,xmax-x],yerr=[(y-ymin),(ymax-y)],fmt='o')
plt.xlabel(r'$E_{eff}(B-V)+const$')
plt.ylabel(r'$\Delta A_B  = A_{eff, B}-3.96 (E_{eff}(B-V)+const)$')
plt.ylim((-0.8,0.1))
pp = PdfPages('output11/Rveffres.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

(y, ymin, ymax) = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k'] ,(50,50-34,50+34),axis=0)
plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
(y, ymin, ymax) = numpy.percentile((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R'] ,(50,50-34,50+34),axis=0)
plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
plt.xlim((-0.1,0.1))
plt.show()
(y, ymin, ymax) = numpy.percentile((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]*fit['k']- ((fit['rho1'][:,1]-fit['rho1'][:,2])[:,None]*fit['R']) ,(50,50-34,50+34),axis=0)
plt.errorbar(x,y,xerr=[x-xmin,xmax-x],yerr=[y-ymin,ymax-y],fmt='o')
plt.show()
x_=numpy.array([-0.08,0.36])
plt.plot(x_,3.96*x_)
plt.xlim((-0.1,0.4))
plt.xlabel(r'$E_{eff}(B-V)+const$')
plt.ylabel(r'$A_{eff, B}$')
pp = PdfPages('output11/ebvebv.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()





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
# pp = PdfPages('output11/Rveff.pdf')
# plt.xlabel(r'$R_{V,eff}$')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# plt.scatter(numpy.median(ab-av,axis=0), numpy.median(rv,axis=0))
# plt.xlabel(r'$E^o(B-V)$')
# plt.ylabel(r'$R_{V,eff}$')
# pp = PdfPages('output11/EVRveff.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
pp = PdfPages('output11/c_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages('output11/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['alpha']-fit['alpha'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$"])
# pp = PdfPages('output11/alpham4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()



figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages('output11/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


figure = corner.corner(fit['eta'],labels=[r"${\eta}_0$",r"${\eta}_1$",r"${\eta}_2$",r"${\eta}_3$",r"${\eta}_4$"])
pp = PdfPages('output11/eta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['beta']-fit['beta'][:,4][:,None]

# figure = corner.corner(mega[:,:4],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$"])
# pp = PdfPages('output11/betam4_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()


figure = corner.corner(fit['gamma'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_2$",r"${\gamma}_3$",r"${\gamma}_4$"])
pp = PdfPages('output11/gamma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['rho1'],labels=[r"${\rho}_{10}$",r"${\rho}_{11}$",r"${\rho}_{12}$",r"${\rho}_{13}$",r"${\rho}_{14}$"])
pp = PdfPages('output11/rho_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2]))[:,None]

# mega = mega-mega[:,2][:,None]
# mega = numpy.delete(mega,2,axis=1)
# mega = numpy.delete(mega,1,axis=1)
# figure = corner.corner(mega,labels=[r"$R_0$",r"$R_3$",r"$R_4$"])
# pp = PdfPages('output11/rx-rv_corner.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages('output11/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

plt.hist(fit['lp__'],cumulative=True) #
plt.xlabel(r'$\ln{p}$')
pp = PdfPages('output11/lp.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


with PdfPages('output11/multipage_pdf.pdf') as pdf:

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
    pp = PdfPages('output11/rx_corner.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()

    mega = fit['rho1']/((fit['rho1'][:,1]-fit['rho1'][:,2]))[:,None]

    figure = corner.corner(mega,labels=[r"$R_{\delta 0}$",r"$R_{\delta 1}$",r"$R_{\delta 2}$",r"$R_{\delta 3}$",r"$R_{\delta 4}$"], \
        range=[[-8.,0.5] for x in xrange(5)])
    pp = PdfPages('output11/rxdelta_corner.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()

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





nlinks = fit['gamma'].shape[0]
mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['eta'],fit['gamma'],fit['rho1'],fit['L_sigma']])
mega = numpy.transpose(mega)

rvs = [1.5,2.5,3.1]
filts = ['U','B','V','R','I']


for index in xrange(5):
    figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index), r"$\alpha_{}$".format(index),\
                    r"$\beta_{}$".format(index),r"$\eta_{}$".format(index),r"$\gamma_{}$".format(index),\
                    r"$\delta_{{1{}}}$".format(index), r"$\sigma_{}$".format(index)])
    figure.suptitle(filts[index],fontsize=28)
    pp = PdfPages('output11/coeff{}.pdf'.format(index))
    plt.savefig(pp,format='pdf')
    pp.close()
    plt.close()



lambdas = numpy.arange(3000.,9000,100)
for rv in rvs:

    f99 = sncosmo.F99Dust(r_v =rv)
    f99.set(ebv=1.37)
    A_ = f99.propagate(lambdas,1.)
    A_=-2.5*numpy.log10(A_)

    norm = f99.propagate([efflam[1]],1.) / f99.propagate(numpy.array([efflam[2]]),1.)
    norm  = -2.5*numpy.log10(norm)
    # A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
    # norm  = sncosmo._extinction.ccm89(numpy.array([efflam[1]]), 1., rv)-sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
    A_ = A_/norm[0]
    plt.plot(lambdas,A_,label=r"$R^F_V={:.1f}$".format(rv))

(y, ymin, ymax) = numpy.percentile(fit['gamma']/((fit['gamma'][:,1]-fit['gamma'][:,2])[:,None]),(50,50-34,50+34),axis=0)
plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
plt.legend()
plt.xlabel(r'Wavelength (\AA)')
plt.ylabel(r'$R_X=\frac{\gamma_X}{\gamma_1-\gamma_2}$')
pp = PdfPages('output11/ccm.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

for rv in rvs:

    f99 = sncosmo.F99Dust(r_v =rv)
    f99.set(ebv=1.37)
    A_ = f99.propagate(lambdas,1.)

    norm  = f99.propagate(numpy.array([efflam[2]]),1.) 
    # A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
    # norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
    # A_ = A_/norm[0]
    A_ = numpy.log10(A_)/numpy.log10(norm)
    plt.plot(lambdas,A_,label=r"$R^F_V={:.1f}$".format(rv))

(y, ymin, ymax) = numpy.percentile(fit['gamma']/fit['gamma'][:,2][:,None],(50,50-34,50+34),axis=0)
plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
plt.legend()
plt.xlabel(r'Wavelength (\AA)')
plt.ylabel(r'$\frac{\gamma_X}{\gamma_2}$')
pp = PdfPages('output11/ccm2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

mega = numpy.array([fit['Delta'].flatten(),fit['EW'][:,:,0].flatten(),fit['EW'][:,:,1].flatten(),fit['sivel'].flatten(), \
    ((fit['gamma'][:,1] - fit['gamma'][:,2])[:,None]*fit['k']).flatten(),((fit['rho1'][:,1] - fit['rho1'][:,2])[:,None]*fit['R']).flatten()])

mega = numpy.transpose(mega)

figure = corner.corner(mega,labels=[r"$\Delta$",r"$EW_{Ca}-109$\AA",r"$EW_{Si}-14$\AA",r"$\lambda_{Si}$",r"$E_\gamma(B-V)$",r"$E_\delta(B-V)$"],range=numpy.zeros(6)+1.)
pp = PdfPages('output11/perobject_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# fig, axes = plt.subplots(nrows=len(filts),sharex=True)
# xerr = numpy.sqrt(EW_cov[:,0,0]) 
# for i in xrange(len(filts)):
#     r = numpy.array( [numpy.argmin(EW_obs[:,0]), numpy.argmax(EW_obs[:,0])])
#     yerr= numpy.sqrt(mag_cov[:,i,i])
#     axes[i].errorbar(EW_obs[:,0],mag_obs[:, i], \
#         xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
#     axes[i].set_ylabel(r'{}'.format(filts[i]))
#     offset  = numpy.mean(mag_obs[:, i] - numpy.median(fit['alpha'][:,i][:,None]*fit['EW'][:,:,0],axis=0))
#     axes[i].plot(EW_obs[r,0],offset+numpy.median(fit['alpha'][:,i],axis=0)*EW_renorm[r,0] \
#         ,color='red',linewidth=2)
# axes[len(filts)-1].set_xlabel(r'EW(Ca)')
# fig.subplots_adjust(hspace=0.001)
# pp = PdfPages('output11/speccamag.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# fig, axes = plt.subplots(nrows=len(filts),sharex=True)
# xerr = numpy.sqrt(EW_cov[:,1,1])
# for i in xrange(len(filts)):
#     r = numpy.array( [numpy.argmin(EW_obs[:,1]), numpy.argmax(EW_obs[:,1])])
#     yerr= numpy.sqrt(mag_cov[:,i,i])
#     axes[i].errorbar(EW_obs[:,1],mag_obs[:, i], \
#         xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
#     axes[i].set_ylabel(r'{}'.format(filts[i]))
#     offset  = numpy.mean(mag_obs[:, i] - numpy.median(fit['beta'][:,i][:,None]*fit['EW'][:,:,1],axis=0))
#     axes[i].plot(EW_obs[r,1],offset+numpy.median(fit['beta'][:,i],axis=0)*EW_renorm[r,1] \
#         ,color='red',linewidth=2)
# axes[len(filts)-1].set_xlabel(r'EW(Si)')
# fig.subplots_adjust(hspace=0.001)
# pp = PdfPages('output11/specsimag.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# fig, axes = plt.subplots(nrows=len(filts),sharex=True)
# xerr = numpy.sqrt(mag_cov[:,1,1]+ mag_cov[:,2,2] - 2*mag_cov[:,1,2])
# for i in xrange(len(filts)):
#     r = numpy.array( [numpy.argmin(mag_obs[:,1]-mag_obs[:,2]), numpy.argmax(mag_obs[:,1]-mag_obs[:,2])])
#     yerr= numpy.sqrt(mag_cov[:,i,i])
#     axes[i].errorbar(mag_obs[:,1]-mag_obs[:,2],mag_obs[:, i], \
#         xerr=[xerr,xerr], yerr=[yerr,yerr],fmt='.')
#     axes[i].set_ylabel(r'{}'.format(filts[i]))
#     slope = numpy.median(fit['gamma'][:,i] / (fit['gamma'][:,1]-fit['gamma'][:,2]),axis=0)
#     offset  = numpy.mean(mag_obs[:, i] - numpy.median(slope*(mag_obs[:,1]-mag_obs[:,2])))
#     axes[i].plot(mag_obs[r,1]-mag_obs[r,2],offset+slope*(mag_obs[r,1]-mag_obs[r,2]) \
#         ,color='red',linewidth=2)
# axes[len(filts)-1].set_xlabel(r'B-V')
# fig.subplots_adjust(hspace=0.001)
# pp = PdfPages('output11/colormag.pdf')
# plt.savefig(pp,format='pdf')
# pp.close()
# plt.close()

# for index in xrange(5):
#     dum = numpy.swapaxes(numpy.array([fit['alpha'][:,index]*color_std,fit['beta'][:,index]*color_std]),0,1)
#     figure = corner.corner(dum,labels=[r"$\alpha_{}$".format(index),r"$\beta_{}$".format(index)],truths=[0,0])
#     pdf.savefig()
#     plt.close()

# plt.plot(fit['ebeta_inv'])
# plt.title('ebeta_inv')
# pdf.savefig()
# plt.close()