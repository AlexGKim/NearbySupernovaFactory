import pickle
import numpy
import pystan


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)
pkl_file.close()

#  The ordering is 'Ca','Si','U','B','V','R','I'

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]
EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

nsne, nmags = mag_obs.shape
# color_obs = numpy.zeros((nsne,nmags))
# color_obs[:,0] = mag_obs[:,0]- mag_obs[:,2]
# color_obs[:,1] = mag_obs[:,1]- mag_obs[:,2]
# color_obs[:,2] = mag_obs[:,2]- mag_obs[:,3]
# color_obs[:,3] = mag_obs[:,2]- mag_obs[:,4]
# color_obs[:,4] = mag_obs[:,2]



# trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1],[0.,0,1,0,0]]
# trans = numpy.array(trans)
# color_cov = numpy.zeros((nsne,5,5))
# for i in xrange(nsne):
#     color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))


# # renormalize the data
EW_mn = EW_obs.mean(axis=0)
# EW_std = EW_obs.std(axis=0)
# EW_renorm = (EW_obs - EW_mn)/EW_std
EW_renorm = (EW_obs - EW_mn)

# # EW_cov_renorm = numpy.array(EW_cov)
# # for i in xrange(nsne):
# #     EW_cov_renorm[i] = EW_cov_renorm[i] /EW_std[None,:] / EW_std[:,None]

# color_mn = color_obs.mean(axis=0)
# # color_std = color_obs[1].std()

# # color_renorm = (color_obs - color_mn)/color_std
# color_renorm = (color_obs - color_mn)
# # color_cov_renorm = numpy.array(color_cov)
# # for i in xrange(nsne):
# #     color_cov_renorm[i] = color_cov_renorm[i] /color_std**2

mag_mn = mag_obs.mean(axis=0)
mag_renorm  = mag_obs-mag_mn

# data = {'D': nsne, 'N_colors': 5, 'N_EWs': 2, 'color_obs': color_renorm, \
#     'EW_obs': EW_renorm, 'EW_cov': EW_cov_renorm, 'color_cov':color_cov_renorm, 'color_std':color_std.astype('float'), 'color_mn':color_mn}
data = {'D': nsne, 'N_mags': 5, 'N_EWs': 2, 'mag_obs': mag_renorm, 'EW_obs': EW_renorm, 'EW_cov': EW_cov, 'mag_cov':mag_cov}

k_simplex = numpy.random.random(size = nsne)
k_simplex /= sum(k_simplex)

init = [{'EW' : EW_renorm, \
    'r_c': numpy.random.uniform(numpy.sqrt(5)-0.01,numpy.sqrt(5)+0.01),\
    'phi_c': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1, size= 3),\
    'phi_c_4': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1),\
    'r_alpha': numpy.random.uniform(numpy.sqrt(5)-0.001,numpy.sqrt(5)+0.001),\
    'phi_alpha': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1, size= 3),\
    'phi_alpha_4': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1),\
    'r_beta': numpy.random.uniform(numpy.sqrt(5)+0.03,numpy.sqrt(5)+0.08),\
    'phi_beta': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1, size= 3),\
    'phi_beta_4': numpy.random.uniform(numpy.pi/4-0.1, numpy.pi/4+0.1),\
    'gamma':numpy.array([1.4,1.21,0.824,0.67])+ numpy.random.uniform(-0.01,0.01,size=4), \
    'k_scale':0.1, 'k_simplex': k_simplex, 'mag_int': mag_renorm+numpy.random.uniform(-0.02,0.02,size=(nsne,5)), \
    'L_sigma':numpy.random.uniform(0.04, 0.08,size=5),'L_Omega':numpy.identity(5)} \
    for _ in range(4)]

# init1 = {'EW' : EW_renorm, 'c': numpy.zeros(5), 'alpha': numpy.zeros(5), \
# 'beta':numpy.array([0.06,0.05,0.04,0.03,0.02]),\
# 'gamma':numpy.array([1.4,1.2,0.85,0.7]), 'k':numpy.zeros(nsne-1), 'mag_int':\
# mag_renorm, 'L_sigma_color':numpy.zeros(5)+0.1,'L_Omega':numpy.identity(5)}

sm = pystan.StanModel(file='gerard3.stan')
# fit = sm.sampling(data=data, iter=200, chains=4,init=[init1,init1,init1,init1])
# control = {'stepsize':1.5}
fit = sm.sampling(data=data, iter=200, chains=4,init=init)
print fit

output = open('temp3.pkl','wb')
pickle.dump(fit.extract(), output)
output.close()
