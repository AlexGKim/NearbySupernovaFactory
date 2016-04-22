import pickle
import numpy
import pystan


pkl_file = open('gege_data.pkl', 'r')
data = pickle.load(pkl_file)

#  The ordering is 'Ca','Si','U','B','V','R','I' 

EW_obs = data['obs'][:,0:2]
mag_obs = data['obs'][:,2:]

nsne, nmags = mag_obs.shape
color_obs = numpy.zeros((nsne,nmags-1))
color_obs[:,0] = mag_obs[:,0]- mag_obs[:,2]
color_obs[:,1] = mag_obs[:,1]- mag_obs[:,2]
color_obs[:,2] = mag_obs[:,2]- mag_obs[:,3]
color_obs[:,3] = mag_obs[:,2]- mag_obs[:,4]

EW_cov = data['cov'][:,0:2,0:2]
mag_cov = data['cov'][:,2:,2:]

pkl_file.close()

trans = [[1.,0,-1,0,0],[0.,1,-1,0,0],[0.,0,1,-1,0],[0.,0,1,0,-1]]
trans = numpy.array(trans)
color_cov = numpy.zeros((nsne,4,4))
for i in xrange(nsne):
    color_cov[i] = numpy.dot(trans,numpy.dot(mag_cov[i], trans.T))


# renormalize the data
EW_mn = EW_obs.mean(axis=0)
EW_std = EW_obs.std(axis=0)

EW_renorm = (EW_obs - EW_mn)/EW_std
EW_cov_renorm = numpy.array(EW_cov)
for i in xrange(nsne):
    EW_cov_renorm[i] = EW_cov_renorm[i] /EW_std[None,:] / EW_std[:,None]

color_mn = color_obs.mean(axis=0)
color_std = color_obs.std(axis=0)

color_renorm = (color_obs - color_mn)/color_std
color_cov_renorm = numpy.array(color_cov)
for i in xrange(nsne):
    color_cov_renorm[i] = color_cov_renorm[i] /color_std[None,:] / color_std[:,None]


data = {'D': nsne, 'N_colors': 4, 'N_EWs': 2, 'color_obs': color_renorm, \
    'EW_obs': EW_renorm, 'EW_cov': EW_cov_renorm, 'color_cov':color_cov_renorm}

sm = pystan.StanModel(file='gerard.stan')
fit = sm.sampling(data=data, iter=200, chains=4)

output = open('temp.pkl','wb')
pickle.dump(fit.extract(), output)
output.close()
