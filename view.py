import pickle
import pystan
import matplotlib.pyplot as plt
from matplotlib import rc
import corner
from matplotlib.backends.backend_pdf import PdfPages
import numpy

f = open('temp.pkl','rb')
fit = pickle.load(f)


with PdfPages('multipage_pdf.pdf') as pdf:

    # plt.plot(fit['EW0'])
    # plt.title('EW0')
    # pdf.savefig()
    # plt.close()

    # plt.plot(fit['L_sigma_EW'])
    # plt.title('L_sigma_EW')
    # pdf.savefig()
    # plt.close()
    
    plt.plot(fit['EW'][:,:10,0])
    plt.title('EW - 0')
    pdf.savefig()
    plt.close()
    
    plt.plot(fit['EW'][:,:10,1])
    plt.title('EW - 1')
    pdf.savefig()
    plt.close()

    plt.plot(fit['c_unnorm'])
    plt.title('c')
    pdf.savefig()
    plt.close()
    
    plt.plot(fit['alpha'])
    plt.title('alpha')
    pdf.savefig()
    plt.close()
    
    plt.plot(fit['beta'])
    plt.title('beta')
    pdf.savefig()
    plt.close()
    
    plt.plot(fit['gamma']*3.1)
    plt.title('gamma')
    pdf.savefig()
    plt.close()
    

    # plt.plot(fit['L_sigma_color'])
    # plt.title('L_sigma_color')
    # pdf.savefig()
    # plt.close()
    
    plt.plot(fit['k'][:,:10]/3.1)
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
    edges = numpy.array([[3300, 4102], [4102, 5100], [5200, 6289], [6289, 7607], [7607, 9200]])
    efflam = []
    for edge in edges:
        efflam.append((edge[1]-edge[0])/2+edge[0])
# [3701, 4601, 5744, 6948, 8403]
    rc('text', usetex=True)
    dif = numpy.array([1.545252680138557- 0.9512179494794494, 1.251613863650247- 0.9512179494794494, 0.9512179494794494-0.7582486105024363, 0.9512179494794494-0.5429015873250165]) 
    dif = dif/dif[1]*3.1
    dum = [dif[0], dif[2], dif[3]]
    figure = corner.corner(fit['gamma']*3.1,labels=[r"${\gamma}_0$",r"${\gamma}_2$",r"${\gamma}_3$"],truths=dum)
    pdf.savefig()
    plt.close()

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