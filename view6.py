import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo

f = open('temp6.pkl','rb')
fit = pickle.load(f)

for key in fit.keys():
    print key, fit[key].min(), fit[key].max()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()


EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]



med = numpy.median(mag_obs[:,1]-mag_obs[:,2])
red = (mag_obs[:,1]-mag_obs[:,2]) > med

EW_obs = EW_obs[red,:]
mag_obs = mag_obs[red,:]
EW_cov = EW_cov[red,:,:]
mag_cov = mag_cov[red,:,:]


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

correction = [ fit['c'][:,i][:,None] + fit['alpha'][:,i][:,None]*fit['EW'][:,:, 0] \
    + fit['beta'][:,i][:,None]*fit['EW'][:,:, 1] \
    + fit['gamma'][:,i][:,None]*fit['k'] \
    + fit['Delta'] \
    for i in xrange(5)]

correction = numpy.array(correction)
correction_median = numpy.median(correction,axis=1)

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
pp = PdfPages('output6/Delta_hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
pp = PdfPages('output6/alpha_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
pp = PdfPages('output6/beta_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['gamma_'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_3$",r"${\gamma}_4$"])
pp = PdfPages('output6/gamma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

# mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages('output6/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
pp = PdfPages('output5/c_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

with PdfPages('output6/multipage_pdf.pdf') as pdf:

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
    # plt.title(r'$\Delta scale$')
    # pdf.savefig()
    # plt.close()

    lineobjects = plt.plot(fit['c'][::10,:],label=['U','B','V','R','I'])
    plt.title('c')
    plt.legend(iter(lineobjects),('U','B','V','R','I'))
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
 
    # plt.plot(fit['R'][::10,:10])
    # plt.title('R')
    # pdf.savefig()
    # plt.close()

    mega = numpy.array([fit['c'],fit['alpha'],fit['beta'],fit['L_sigma']])
    mega = numpy.transpose(mega)

    for index in xrange(5):
        figure = corner.corner(mega[index,:,:],labels=[r"$c_{}$".format(index),r"$\alpha_{}$".format(index),\
            r"$\beta_{}$".format(index),r"$\sigma_{}$".format(index)])
        pdf.savefig()
        plt.close()

    # filter wavelengths
edges = numpy.array([[3300., 4102], [4102, 5100], [5200, 6289], [6289, 7607], [7607, 9200]])
efflam = []
for edge in edges:
    efflam.append((edge[1]-edge[0])/2+edge[0])
# [3701, 4601, 5744, 6948, 8403]
rc('text', usetex=True)


rvs = [1.,3.1,4.1]
filts = ['U','B','V','R','I']

fig, axes = plt.subplots(nrows=2, ncols=len(filts))
for i in xrange(len(filts)):

    axes[0,i].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, i])
    axes[0,i].set_xlabel(r'B-V')
    axes[0,i].set_ylabel(r'{}'.format(filts[i]))

    axes[1,i].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, i]-correction_median[i,:])
    axes[1,i].set_xlabel(r'B-V')
    axes[1,i].set_ylabel(r'{} + correction'.format(filts[i]))


plt.tight_layout()
pp = PdfPages('output6/colormag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

fig, axes = plt.subplots(nrows=2, ncols=len(filts))
for i in xrange(len(filts)):
    axes[0,i].scatter(EW_obs[:,1],mag_obs[:, i])
    axes[0,i].set_xlabel(r'EW(Si)')
    axes[0,i].set_ylabel(r'{}'.format(filts[i]))

    axes[1,i].scatter(EW_obs[:,1],mag_obs[:, i]-correction_median[i,:])
    axes[1,i].set_xlabel(r'EW(Si)')
    axes[1,i].set_ylabel(r'{} + correction'.format(filts[i]))

plt.tight_layout()
pp = PdfPages('output6/specsimag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

fig, axes = plt.subplots(nrows=2, ncols=len(filts))
for i in xrange(len(filts)):


    axes[0,i].scatter(EW_obs[:,0],mag_obs[:, i])
    axes[0,i].set_xlabel(r'EW(Ca)')
    axes[0,i].set_ylabel(r'{}'.format(filts[i]))

    axes[1,i].scatter(EW_obs[:,0],mag_obs[:, i]-correction_median[i,:])
    axes[1,i].set_xlabel(r'EW(Ca)')
    axes[1,i].set_ylabel(r'{} + correction'.format(filts[i]))



plt.tight_layout()
pp = PdfPages('output6/speccamag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

lambdas = numpy.arange(3500.,9000,100)
for rv in rvs:
    A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
    norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
    A_ = A_/norm[0]
    plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))

(y, ymin, ymax) = numpy.percentile(fit['gamma'],(50,32,64),axis=0)plt.errorbar(efflam,y,yerr=[y-ymin,ymax-y],fmt='o')
plt.legend()
plt.xlabel(r'Wavelength (\AA)')
pp = PdfPages('output6/ccm.pdf')
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