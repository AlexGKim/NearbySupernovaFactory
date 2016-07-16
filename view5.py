import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo

color_std=0.272560118531

f = open('temp5.pkl','rb')
fit = pickle.load(f)

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

# dum = numpy.zeros((5,5))
# for x1, x2 in zip(fit['L_Omega'], fit['L_sigma']):
#     dum= dum+ numpy.dot(x2[:,None],x2[None,:])*numpy.dot(x1,x1.T)

# dum/= len(fit['L_Omega'])
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in dum])

# trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
# trans = numpy.array(trans)
# color_cov = numpy.dot(trans,numpy.dot(dum, trans.T))
# print " \\\\\n".join([" & ".join(map('{0:.4f}'.format, line)) for line in color_cov])



gamma_median = numpy.median(fit['gamma'],axis=0)
rho0_median = numpy.median(fit['rho0'],axis=0)
rho1_median = numpy.median(fit['rho1'],axis=0)

R_median = numpy.median(fit['R'],axis=0)
EW_median = numpy.median(fit['EW'],axis=0)
beta_median = numpy.median(fit['beta'],axis=0)
alpha_median = numpy.median(fit['alpha'],axis=0)
k_median =  numpy.median(fit['k'],axis=0)
Delta_median =  numpy.median(fit['Delta'],axis=0)
R_median =  numpy.median(fit['R'],axis=0)

# filtname = ['U','B','V','R','I']

# for i, fi in enumerate(filtname):
#     plt.hist([mag_obs[:, i], mag_obs[:, i] - gamma_median[i]*k_median, \
#         mag_obs[:, i] - gamma_median[i]*k_median -beta_median[i]*EW_median[:,1], \
#         mag_obs[:, i] - gamma_median[i]*k_median -beta_median[i]*EW_median[:,1]-alpha_median[i]*EW_median[:,0]], \
#         label=[r'$M_o$', r'$M_o - \gamma k$', r'$M_o - \gamma k - \beta EW_{Si}$', \
#         r'$M_o - \gamma k - \beta EW_{Si} - \alpha EW_{Ca}$'])
#     plt.legend()
#     plt.title('Filter {}'.format(fi))
#     pp = PdfPages('output5/magresidual_{}.pdf'.format(fi))
#     plt.savefig(pp, format='pdf')
#     pp.close()
#     plt.close()

# temp = (fit['Delta']-1./fit['Delta'].shape[1])*fit['Delta_scale'][:,None]
# # delta_median =  numpy.median(temp,axis=0)
# print temp.flatten().std()
print fit['Delta'].flatten().std()
plt.hist(fit['Delta'].flatten(),normed=True,bins=20)
plt.title(r'$\Delta$')
pp = PdfPages('output5/Delta_hist.pdf')
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

figure = corner.corner(fit['rho0'],labels=[r"$\rho_{00}$",r"$\rho_{01}$",r"$\rho_{02}$",r"$\rho_{03}$",r"$\rho_{04}$"])
pp = PdfPages('output5/rho0_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

figure = corner.corner(fit['rho1'],labels=[r"$\rho_{10}$",r"$\rho_{11}$",r"$\rho_{12}$",r"$\rho_{13}$",r"$\rho_{13}$"])
pp = PdfPages('output5/rho1_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()


mega = numpy.concatenate((fit['Delta_scale'][:,None],fit['L_sigma']),axis=1)
figure = corner.corner(mega,labels=[r"$\Delta$ scale",r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
pp = PdfPages('output5/sigma_corner.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()

with PdfPages('output5/multipage_pdf.pdf') as pdf:

    # plt.plot(fit['EW0'])
    # plt.title('EW0')
    # pdf.savefig()
    # plt.close()

    # lineobjects = plt.plot(fit['sigma_Delta'][::10],label=[])
    # plt.title(r'$\sigma_\Delta$')
    # pdf.savefig()
    # plt.close()
    lineobjects = plt.plot(fit['lp__'][::10])
    plt.title(r'log p')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(fit['rho0'][::10,:])
    plt.title(r'$\rho_0$ ')
    plt.legend(iter(lineobjects),(r'$\rho_0$',r'$\rho_1$',r'$\rho_3$'))
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(fit['rho1'][::10,:])
    plt.title(r'$\rho_1$')
    plt.legend(iter(lineobjects),(r'$\rho_0$',r'$\rho_1$',r'$\rho_2$',r'$\rho_3$'))
    pdf.savefig()
    plt.close()

    # plt.plot(fit['R_mn'][::10])
    # plt.title(r'R_mn')
    # pdf.savefig()
    # plt.close()

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
    plt.title(r'$\Delta scale$')
    pdf.savefig()
    plt.close()

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
 

    # plt.plot(fit['L_sigma_color'])
    # plt.title('L_sigma_color')
    # pdf.savefig()
    # plt.close()
    
    # plt.plot(fit['k_simplex'][::10,:10])
    # plt.title('k')
    # pdf.savefig()
    # plt.close()
    
    plt.plot(fit['R'][::10,:10])
    plt.title('R')
    pdf.savefig()
    plt.close()
    


    # plt.plot(fit['EXY'][:,:10,1])
    # plt.title('EXY')
    # pdf.savefig()
    # plt.close()
    
    # plt.plot(fit['colors'][:,:10,1])
    # plt.title('colors')
    # pdf.savefig()
    # plt.close()

    # dif = numpy.array([1.545252680138557- 0.9512179494794494, 1.251613863650247- 0.9512179494794494, 0.9512179494794494-0.7582486105024363, 0.9512179494794494-0.5429015873250165,0.9512179494794494])
    # dif = dif/dif[1]

    figure = corner.corner(fit['c'],labels=[r"${c}_0$",r"${c}_1$",r"${c}_2$",r"${c}_3$",r"${c}_4$"])
    pdf.savefig()
    plt.close()



    # plt.hist(fit['k_simplex'].flatten())
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



fig, axes = plt.subplots(nrows=2, ncols=2)
rvs = [1.,3.1,4.1]
#fig.subplots_adjust(hspace=1, vspace=0.5)

axes[0,0].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, 0])
xl = numpy.array(axes[0,0].get_xlim()) * numpy.array([0.9, 0.6])
yl = numpy.array(axes[0,0].get_ylim())
for rv in rvs:
    A_ = sncosmo._extinction.ccm89(numpy.array(efflam), 1., rv)
    slope = A_[0]/(A_[1]-A_[2])
    y0 = numpy.mean(mag_obs[:, 0] - slope*(mag_obs[:, 1]-mag_obs[:, 2]))
    axes[0,0].plot(xl,slope*xl +y0,label=rv)
    axes[0,0].legend(loc=4)
    axes[0,0].set_xlabel(r'B-V')
    axes[0,0].set_ylabel(r'U')

axes[1,0].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, 0]-beta_median[0]*EW_median[:,1]-alpha_median[0]*EW_median[:,0]-rho0_median[0]*R_median-rho1_median[0]*R_median**2)
xl = numpy.array(axes[1,0].get_xlim())* numpy.array([0.9, 0.6])
yl = numpy.array(axes[1,0].get_ylim())
for rv in rvs:
    A_ = sncosmo._extinction.ccm89(numpy.array(efflam), 1., rv)
    slope = A_[0]/(A_[1]-A_[2])
    y0 = numpy.mean(mag_obs[:, 0] -beta_median[0]*EW_median[:,1]-alpha_median[0]*EW_median[:,0] -rho0_median[0]*R_median -rho1_median[0]*R_median**2- slope*(mag_obs[:, 1]-mag_obs[:, 2]))
    axes[1,0].plot(xl,slope*xl +y0,label=rv)
    axes[1,0].legend(loc=4)
    axes[1,0].set_xlabel(r'B-V')
    axes[1,0].set_ylabel(r'U + spec correction')

axes[0,1].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, 4])
xl = numpy.array(axes[0,0].get_xlim())* numpy.array([0.9, 0.6])
yl = numpy.array(axes[0,0].get_ylim())

for rv in rvs:
    A_ = sncosmo._extinction.ccm89(numpy.array(efflam), 1., rv)
    slope = A_[4]/(A_[1]-A_[2])
    y0 = numpy.mean(mag_obs[:, 4] - slope*(mag_obs[:, 1]-mag_obs[:, 2]))
    axes[0,1].plot(xl,slope*xl +y0,label=rv)
    axes[0,1].legend(loc=4)
    axes[0,1].set_xlabel(r'B-V')
    axes[0,1].set_ylabel(r'I')

axes[1,1].scatter(mag_obs[:, 1]-mag_obs[:, 2],mag_obs[:, 4]-beta_median[4]*EW_median[:,1]-alpha_median[4]*EW_median[:,0])
xl = numpy.array(axes[1,1].get_xlim())* numpy.array([0.9, 0.6])
yl = numpy.array(axes[1,1].get_ylim())

for rv in rvs:
    A_ = sncosmo._extinction.ccm89(numpy.array(efflam), 1., rv)
    slope = A_[4]/(A_[1]-A_[2])
    y0 = numpy.mean(mag_obs[:, 4] -beta_median[4]*EW_median[:,1]-alpha_median[4]*EW_median[:,0] - slope*(mag_obs[:, 1]-mag_obs[:, 2]))
    axes[1,1].plot(xl,slope*xl +y0,label=rv)
    axes[1,1].legend(loc=4)
    axes[1,1].set_xlabel(r'B-V')
    axes[1,1].set_ylabel(r'I + spec correction')

plt.tight_layout()
pp = PdfPages('output5/colormag.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.close()



# lambdas = numpy.arange(3500.,9000,100)
# for rv in rvs:
#     A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
#     norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
#     A_ = A_/norm[0]
#     plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))

# dum = (1-.68)/2
# gamma_sort = numpy.sort(fit['gamma'],axis=0)
# ngamma = gamma_sort.shape[0]
# y = gamma_sort[ngamma/2,:]
# ymin = gamma_sort[ngamma*dum,:]
# ymax  = gamma_sort[ngamma*(1-dum),:]
# plt.errorbar(numpy.delete(efflam,2),y,yerr=[y-ymin,ymax-y],fmt='o')
# plt.legend()
# plt.xlabel(r'Wavelength (\AA)')
# pp = PdfPages('output5/ccm.pdf')
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
    


#   A_U 1.545252680138557, 1.25, 0.952, 0.758, 0.543

# cov = numpy.mean(fit['Omega_color'],axis=0)
# print cov

# plt.hist(fit['c'][:,:])


# plt.clf()
# plt.plot(fit['k'][:,0:10])
# plt.hist(fit['EXY'][900,:,1])

# figure = corner.corner(fit['c'])
# figure.savefig("corner.png")
