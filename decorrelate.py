#!/usr/bin/env python

import pickle
import numpy
import copy

def fix3_x1():
    f = open('fix3_x1.pkl','rb')
    (fit,_) = pickle.load(f)

    outfit= copy.deepcopy(fit)

    for i in xrange(fit['x1'].shape[0]):
    # for i in xrange(1):
        # for j in xrange(fit['Delta'].shape[1]):
        #     print "original"
        #     print fit['c'][i,:]+fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * fit['k'][i,j]  + fit['rho1'][i,:]* fit['R'][i,j]

        #     newk = fit['k'][i,j] + fit['EW'][i,j,0]*epsilon
        #     newalpha  = fit['alpha'][i,:]-fit['gamma'][i,:]*epsilon
        #     # newc = fit['c'][i,:] +  fit['EW'][i,:,0].mean()*epsilon
        #     print fit['c'][i,:]+  newalpha *fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * newk  + fit['rho1'][i,:]* fit['R'][i,j]

        hs = numpy.array([fit['EW'][i,:,0],fit['EW'][i,:,1],fit['sivel'][i,:], fit['x1'][i,:], \
            fit['k'][i,:], fit['R'][i,:], fit['ev_sig'][i]*fit['mag_int_raw'][i,:],fit['Delta'][i,:]])
        bigcov = numpy.cov(hs)
        hscov = bigcov[:7,:7]
        hdcov = bigcov[:7,7]
        c = numpy.dot(numpy.linalg.inv(hscov), hdcov)

        newalpha = fit['alpha'][i,:] + c[0]
        newbeta = fit['beta'][i,:] + c[1]
        neweta = fit['eta'][i,:] + c[2]
        newzeta = fit['zeta'][i,:] + c[3]
        newgamma = fit['gamma'][i,:] + c[4]
        newrho1 = fit['rho1'][i,:] + c[5]
        newev = fit['ev'][i,:] + c[6]

        newDelta=[]
        for j in xrange(fit['Delta'].shape[1]):

            # print "original"
            # print fit['Delta'][i,j] + fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
            #     fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j]

            newDelta.append(fit['Delta'][i,j] - \
                numpy.dot(c,numpy.array([fit['EW'][i,j,0],fit['EW'][i,j,1],fit['sivel'][i,j], fit['x1'][i,j], \
                fit['k'][i,j], fit['R'][i,j],fit['ev_sig'][i]* fit['mag_int_raw'][i,j]])))

            # print "new" 
            # print newDelta + newalpha*fit['EW'][i,j,0]+newbeta* fit['EW'][i,j,1] +\
            #     neweta* fit['sivel'][i,j] + newzeta* fit['x1'][i,j]

        newDelta = numpy.array(newDelta)

        # # # now do latent parameters
        # latents =[fit['k'][i,:],fit['R'][i,:],fit['ev_sig'][i]*fit['mag_int_raw'][i,:]]
        # cofactors =[newgamma,newrho1,newev]
        # newpars = []
        # for lat, cof in zip(latents,cofactors):
        #     hs = numpy.array([fit['EW'][i,:,0],fit['EW'][i,:,1],fit['sivel'][i,:], fit['x1'][i,:],lat])
        #     bigcov = numpy.cov(hs)
        #     hscov = bigcov[:4,:4]
        #     hdcov = bigcov[:4,4]
        #     c = numpy.dot(numpy.linalg.inv(hscov), hdcov)

        #     newalpha = newalpha + c[0]*cof
        #     newbeta = newbeta + c[1]*cof
        #     neweta = neweta + c[2]*cof
        #     newzeta = newzeta + c[3]*cof

        #     newk=[]
        #     for j in xrange(fit['Delta'].shape[1]):

        #         # print "original"
        #         # print fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         #     fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] + fit['gamma'][i,:]*fit['k'][i,j]

        #         newk.append(lat[j] - \
        #             numpy.dot(c,numpy.array([fit['EW'][i,j,0],fit['EW'][i,j,1],fit['sivel'][i,j], fit['x1'][i,j]])))
        #         # print newk[-1]
        #         # print "new"
        #         # print newalpha*fit['EW'][i,j,0]+newbeta* fit['EW'][i,j,1] +\
        #         #     neweta* fit['sivel'][i,j] + newzeta* fit['x1'][i,j] + fit['gamma'][i,:]*newk[-1]

        #     newpars.append(newk)

        # newpars = numpy.array(newpars)


        outfit['alpha'][i,:] = newalpha
        outfit['beta'][i,:] = newbeta
        outfit['eta'][i,:]  = neweta
        outfit['zeta'][i,:] = newzeta
        outfit['Delta'][i,:] = newDelta
        outfit['gamma'][i,:] = newgamma
        outfit['rho1'][i,:] = newrho1
        outfit['ev'][i,:] = newev
        # outfit['k'][i,:] = newpars[0]
        # outfit['R'][i,:] = newpars[1]
        # outfit['mag_int_raw'][i,:]  = newpars[2]/fit['ev_sig'][i]

    # for i in xrange(1):#fit['x1'].shape[0]):
    #     for j in xrange(fit['Delta'].shape[1]):
    #         # print ( fit['mag_int_raw'][i,j]," ",outfit['mag_int_raw'][i,j])
    #         term1 = fit['Delta'][i,j] + fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
    #             fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] + fit['gamma'][i,:]*fit['k'][i,j]+\
    #             fit['rho1'][i,:]*fit['R'][i,j] +fit['ev_sig'][i] * fit['ev'][i,:]* fit['mag_int_raw'][i,j]
    #         term2 = outfit['Delta'][i,j] + outfit['alpha'][i,:]*outfit['EW'][i,j,0]+outfit['beta'][i,:]* outfit['EW'][i,j,1] +\
    #             outfit['eta'][i,:]* outfit['sivel'][i,j] + outfit['zeta'][i,:]* outfit['x1'][i,j] + outfit['gamma'][i,:]*outfit['k'][i,j]+\
    #             outfit['rho1'][i,:]*outfit['R'][i,j] +outfit['ev_sig'][i]* outfit['ev'][i,:]* outfit['mag_int_raw'][i,j]
    #         print term1-term2
    #     hs = numpy.array([outfit['EW'][i,:,0],outfit['EW'][i,:,1],outfit['sivel'][i,:], outfit['x1'][i,:],outfit['Delta'][i,:], \
    #         outfit['k'][i,:],outfit['R'][i,:],outfit['mag_int_raw'][i,:]])
    #     print numpy.cov(hs)    

    output = open('fix3_x1_decorr.pkl','wb')
    pickle.dump(outfit, output, protocol=2)
    output.close()


def fix3():
    f = open('fix3.pkl','rb')
    (fit,_) = pickle.load(f)

    outfit= copy.deepcopy(fit)

    for i in xrange(fit['Delta'].shape[0]):
    # for i in xrange(1):
        # for j in xrange(fit['Delta'].shape[1]):
        #     print "original"
        #     print fit['c'][i,:]+fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * fit['k'][i,j]  + fit['rho1'][i,:]* fit['R'][i,j]

        #     newk = fit['k'][i,j] + fit['EW'][i,j,0]*epsilon
        #     newalpha  = fit['alpha'][i,:]-fit['gamma'][i,:]*epsilon
        #     # newc = fit['c'][i,:] +  fit['EW'][i,:,0].mean()*epsilon
        #     print fit['c'][i,:]+  newalpha *fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * newk  + fit['rho1'][i,:]* fit['R'][i,j]

        hs = numpy.array([fit['EW'][i,:,0],fit['EW'][i,:,1],fit['sivel'][i,:],  \
            fit['k'][i,:], fit['R'][i,:], fit['ev_sig'][i]*fit['mag_int_raw'][i,:],fit['Delta'][i,:]])
        bigcov = numpy.cov(hs)
        hscov = bigcov[:6,:6]
        hdcov = bigcov[:6,6]
        c = numpy.dot(numpy.linalg.inv(hscov), hdcov)

        newalpha = fit['alpha'][i,:] + c[0]
        newbeta = fit['beta'][i,:] + c[1]
        neweta = fit['eta'][i,:] + c[2]
        newgamma = fit['gamma'][i,:] + c[3]
        newrho1 = fit['rho1'][i,:] + c[4]
        newev = fit['ev'][i,:] + c[5]

        newDelta=[]
        for j in xrange(fit['Delta'].shape[1]):
            newDelta.append(fit['Delta'][i,j] - \
                numpy.dot(c,numpy.array([fit['EW'][i,j,0],fit['EW'][i,j,1],fit['sivel'][i,j],  \
                fit['k'][i,j], fit['R'][i,j],fit['ev_sig'][i]* fit['mag_int_raw'][i,j]])))

        newDelta = numpy.array(newDelta)

        outfit['alpha'][i,:] = newalpha
        outfit['beta'][i,:] = newbeta
        outfit['eta'][i,:]  = neweta
        outfit['Delta'][i,:] = newDelta
        outfit['gamma'][i,:] = newgamma
        outfit['rho1'][i,:] = newrho1
        outfit['ev'][i,:] = newev

    output = open('fix3_decorr.pkl','wb')
    pickle.dump(outfit, output, protocol=2)
    output.close()


def fix1():
    f = open('fix1.pkl','rb')
    (fit,_) = pickle.load(f)

    outfit= copy.deepcopy(fit)

    for i in xrange(fit['Delta'].shape[0]):
    # for i in xrange(1):
        # for j in xrange(fit['Delta'].shape[1]):
        #     print "original"
        #     print fit['c'][i,:]+fit['alpha'][i,:]*fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * fit['k'][i,j]  + fit['rho1'][i,:]* fit['R'][i,j]

        #     newk = fit['k'][i,j] + fit['EW'][i,j,0]*epsilon
        #     newalpha  = fit['alpha'][i,:]-fit['gamma'][i,:]*epsilon
        #     # newc = fit['c'][i,:] +  fit['EW'][i,:,0].mean()*epsilon
        #     print fit['c'][i,:]+  newalpha *fit['EW'][i,j,0]+fit['beta'][i,:]* fit['EW'][i,j,1] +\
        #         fit['eta'][i,:]* fit['sivel'][i,j] + fit['zeta'][i,:]* fit['x1'][i,j] +\
        #         fit['gamma'][i,:] * newk  + fit['rho1'][i,:]* fit['R'][i,j]

        hs = numpy.array([fit['EW'][i,:,0],fit['EW'][i,:,1],fit['sivel'][i,:],  \
            fit['k'][i,:], fit['R'][i,:],fit['Delta'][i,:]])
        bigcov = numpy.cov(hs)
        hscov = bigcov[:5,:5]
        hdcov = bigcov[:5,5]
        c = numpy.dot(numpy.linalg.inv(hscov), hdcov)

        newalpha = fit['alpha'][i,:] + c[0]
        newbeta = fit['beta'][i,:] + c[1]
        neweta = fit['eta'][i,:] + c[2]
        newgamma = fit['gamma'][i,:] + c[3]
        newrho1 = fit['rho1'][i,:] + c[4]

        newDelta=[]
        for j in xrange(fit['Delta'].shape[1]):
            newDelta.append(fit['Delta'][i,j] - \
                numpy.dot(c,numpy.array([fit['EW'][i,j,0],fit['EW'][i,j,1],fit['sivel'][i,j],  \
                fit['k'][i,j], fit['R'][i,j]])))

        newDelta = numpy.array(newDelta)

        outfit['alpha'][i,:] = newalpha
        outfit['beta'][i,:] = newbeta
        outfit['eta'][i,:]  = neweta
        outfit['Delta'][i,:] = newDelta
        outfit['gamma'][i,:] = newgamma
        outfit['rho1'][i,:] = newrho1

    output = open('fix1_decorr.pkl','wb')
    pickle.dump(outfit, output, protocol=2)
    output.close()

fix1()