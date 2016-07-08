import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo

color_std=0.272560118531

f = open('temp3.pkl','rb')
fit = pickle.load(f)
f.close()

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

with PdfPages('multipage_pdf.pdf') as pdf:

    gamma_median = numpy.median(fit['gamma'],axis=0)
    gamma_median = numpy.insert(gamma_median,2,1.)

    k_median =  numpy.append(0.,numpy.median(fit['k'],axis=0))
    EW_median = numpy.median(fit['EW'],axis=0)
    beta_median = numpy.median(fit['beta'],axis=0)
    alpha_median = numpy.median(fit['alpha'],axis=0)

    filtname = ['U','B','V','R','I']
    for i, fi in enumerate(filtname):
        plt.hist([mag_obs[:, i], mag_obs[:, i] - gamma_median[i]*k_median, \
            mag_obs[:, i] - gamma_median[i]*k_median -beta_median[i]*EW_median[:,1], \
            mag_obs[:, i] - gamma_median[i]*k_median -beta_median[i]*EW_median[:,1]-alpha_median[i]*EW_median[:,0]], \
            label=[r'$M_o$', r'$M_o - \gamma k$', r'$M_o - \gamma k - \beta EW_{Si}$', \
            r'$M_o - \gamma k - \beta EW_{Si} - \alpha EW_{Ca}$'])
        plt.legend()
        plt.title('Filter {}'.format(fi))
        pdf.savefig()
        plt.close()

    # plt.plot(fit['EW0'])
    # plt.title('EW0')
    # pdf.savefig()
    # plt.close()

    lineobjects = plt.plot(fit['lp__'][::10])
    plt.title(r'log p')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(fit['L_sigma'][::10,:])
    plt.title(r'$\sigma$')
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
    
    lineobjects = plt.plot(fit['gamma'][::10],label=['U','B','R','I'])
    plt.title('gamma')
    plt.legend(iter(lineobjects),('U','B','R','I'))
    pdf.savefig()
    plt.close()
    

    # plt.plot(fit['L_sigma_color'])
    # plt.title('L_sigma_color')
    # pdf.savefig()
    # plt.close()
    
    plt.plot(fit['k'][::10,:10])
    plt.title('k')
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

    figure = corner.corner(fit['alpha'],labels=[r"${\alpha}_0$",r"${\alpha}_1$",r"${\alpha}_2$",r"${\alpha}_3$",r"${\alpha}_4$"])
    pdf.savefig()
    plt.close()

    figure = corner.corner(fit['beta'],labels=[r"${\beta}_0$",r"${\beta}_1$",r"${\beta}_2$",r"${\beta}_3$",r"${\beta}_4$"])
    pdf.savefig()
    plt.close()

    figure = corner.corner(fit['gamma'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_3$",r"${\gamma}_4$"])
    pdf.savefig()
    plt.close()

    figure = corner.corner(fit['L_sigma'],labels=[r"${\sigma}_0$",r"${\sigma}_1$",r"${\sigma}_2$",r"${\sigma}_3$",r"${\sigma}_4$"])
    pdf.savefig()
    plt.close()

    dum = numpy.median(fit['k'],axis=0)
    plt.hist(dum)
    plt.title('k')
    pdf.savefig()
    plt.close()

    c_con = numpy.concatenate((fit['r_c'][:,None],fit['phi_c'],fit['phi_c_4'][:,None]),axis=1)
    figure = corner.corner(c_con,labels=["r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"])
    plt.title('c')
    pdf.savefig()
    plt.close()

    plt.plot(c_con[::10,0])
    plt.title(r'$r_c$')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(c_con[::10,:])
    plt.title(r'$c$')
    plt.legend(iter(lineobjects),("r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"))
    pdf.savefig()
    plt.close()

    alpha_con = numpy.concatenate((fit['r_alpha'][:,None],fit['phi_alpha'],fit['phi_alpha_4'][:,None]),axis=1)
    figure = corner.corner(alpha_con,labels=["r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"])
    plt.title('alpha')
    pdf.savefig()
    plt.close()

    plt.plot(alpha_con[::10,0])
    plt.title(r'$r_\alpha$')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(alpha_con[::10,:])
    plt.title(r'$\alpha$')
    plt.legend(iter(lineobjects),("r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"))
    pdf.savefig()
    plt.close()

    beta_con= numpy.concatenate((fit['r_beta'][:,None],fit['phi_beta'],fit['phi_beta_4'][:,None]),axis=1)
    figure = corner.corner(beta_con,labels=["r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"])
    plt.title('beta')
    pdf.savefig()
    plt.close()

    plt.plot(beta_con[::10,0])
    plt.title(r'$r_\beta$')
    pdf.savefig()
    plt.close()

    lineobjects = plt.plot(beta_con[::10,:])
    plt.title(r'$\beta$')
    plt.legend(iter(lineobjects),("r", r"$\phi_0$",r"$\phi_1$",r"$\phi_2$",r"$\phi_3$"))
    pdf.savefig()
    plt.close()


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

    lambdas = numpy.arange(3500.,9000,100)
    rvs = [1.8,3.1,4.8]
    for rv in rvs:
        A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
        norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
        A_ = A_/norm[0]
        plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))

    dum = (1-.68)/2

    ymin = numpy.percentile(fit['gamma'],dum*100,axis=0)
    ymax = numpy.percentile(fit['gamma'],(1-dum)*100,axis=0)
    gamma_median = numpy.median(fit['gamma'],axis=0)
    # gamma_sort[int(round(ngamma*dum)),:]
    # ymax  = gamma_sort[int(round(ngamma*(1-dum))),:]
    plt.errorbar(numpy.delete(efflam,2),gamma_median,yerr=[gamma_median-ymin,ymax-gamma_median],fmt='o')
    plt.legend()
    plt.xlabel(r'wavelength (\AA)')
    plt.ylabel(r'$A(\lambda)/A_V$')
    pdf.savefig()
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
    


#   A_U 1.545252680138557, 1.25, 0.952, 0.758, 0.543

# cov = numpy.mean(fit['Omega_color'],axis=0)
# print cov

# plt.hist(fit['c'][:,:])


# plt.clf()
# plt.plot(fit['k'][:,0:10])
# plt.hist(fit['EXY'][900,:,1])

# figure = corner.corner(fit['c'])
# figure.savefig("corner.png")