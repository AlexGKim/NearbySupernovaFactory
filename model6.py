import pickle
import numpy
import pystan

# 1-color model with bluest half

pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

med = numpy.percentile(mag_obs[:,1]-mag_obs[:,2],50)
blue = (mag_obs[:,1]-mag_obs[:,2]) < med
 
EW_obs = EW_obs[blue,:]
mag_obs = mag_obs[blue,:]
EW_cov = EW_cov[blue,:,:]
mag_cov = mag_cov[blue,:,:]

nsne, nmags = mag_obs.shape

# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_renorm = (EW_obs - EW_mn)

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov}

Delta_simplex = numpy.zeros(nsne)+1./nsne
k_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)
k_simplex = k_simplex/k_simplex.sum()

init = [{'EW' : EW_renorm, \
         'c': numpy.zeros(5),\
         'alpha': numpy.zeros(5), \
         'beta':numpy.zeros(5), \
         'gamma0': 0.1,\
         'gamma_': numpy.zeros(4), \
         'mag_int': mag_renorm, \
         'L_sigma': numpy.zeros(5)+0.05, \
         'L_Omega': numpy.identity(5), \
         'Delta_unit':Delta_simplex, 'Delta_scale': nsne/8.,\
         'k_unit': k_simplex} \
        for _ in range(4)]

sm = pystan.StanModel(file='gerard7.stan')
control = {'stepsize':1.}
fit = sm.sampling(data=data, iter=2000, chains=4,control=control,init=init,thin=2)
print fit

output = open('temp6.pkl','wb')
pickle.dump(fit.extract(), output)
output.close()