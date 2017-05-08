import numpy

def flip(fit):
	wmn, wless, wmore = numpy.percentile(fit['rho1'][:,0][:,None]*fit['R'],(50,50-34,50+34),axis=0)
	shit = wmn/(wmn-wless)
	shit[wmn < 0] = wmn[wmn < 0]/(wmn[wmn < 0]-wmore[wmn < 0])
	wmax = numpy.argmax(shit)
	temprho1 = numpy.array(fit['rho1'])
	tempR = numpy.array(fit['R'])
	for i in xrange(temprho1.shape[0]):
	    if (fit['R'][i,wmax] < 0):
	        temprho1[i,:] *=-1
	        tempR[i,:] *=-1
	return temprho1,tempR