import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sncosmo

color_std=0.272560118531

f = open('temp2.pkl','rb')
fit = pickle.load(f)


with PdfPages('multipage_pdf.pdf') as pdf:

    # plt.plot(fit['EW0'])
    # plt.title('EW0')
    # pdf.savefig()
    # plt.close()

    lineobjects = plt.plot(fit['L_sigma'],label=[])
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

    lineobjects = plt.plot(fit['c'],label=['U','B','V','R','I'])
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
    
    lineobjects = plt.plot(fit['gamma'],label=['U','B','R','I'])
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

    # filter wavelengths
    edges = numpy.array([[3300., 4102], [4102, 5100], [5200, 6289], [6289, 7607], [7607, 9200]])
    efflam = []
    for edge in edges:
        efflam.append((edge[1]-edge[0])/2+edge[0])
# [3701, 4601, 5744, 6948, 8403]
    rc('text', usetex=True)

    lambdas = numpy.arange(3500.,9000,100)
    rvs = [2.1,3.1,4.1]
    for rv in rvs:
        A_ = sncosmo._extinction.ccm89(lambdas, 1., rv)
        norm  = sncosmo._extinction.ccm89(numpy.array([efflam[2]]), 1., rv)
        A_ = A_/norm[0]
        plt.plot(lambdas,A_,label=r"$R_V={:.1f}$".format(rv))

    dum = (1-.68)/2
    gamma_sort = numpy.sort(fit['gamma'],axis=0)
    ngamma = gamma_sort.shape[0]
    y = gamma_sort[ngamma/2,:]
    ymin = gamma_sort[ngamma*dum,:]
    ymax  = gamma_sort[ngamma*(1-dum),:]
    plt.errorbar(numpy.delete(efflam,2),y,yerr=[y-ymin,ymax-y],fmt='o')
    plt.legend()
    pdf.savefig()
    plt.close()
    # dif = numpy.array([1.545252680138557- 0.9512179494794494, 1.251613863650247- 0.9512179494794494, 0.9512179494794494-0.7582486105024363, 0.9512179494794494-0.5429015873250165,0.9512179494794494])
    # dif = dif/dif[1]

    # dum = [dif[0], dif[2], dif[3], dif[4]]
    figure = corner.corner(fit['gamma'],labels=[r"${\gamma}_0$",r"${\gamma}_1$",r"${\gamma}_3$",r"${\gamma}_4$"])
    pdf.savefig()
    plt.close()

    figure = corner.corner(fit['L_sigma'],labels=[r"${L_\sigma}_0$",r"${L_\sigma}_1$",r"${L_\sigma}_2$",r"${L_\sigma}_3$",r"${L_\sigma}_4$"])
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