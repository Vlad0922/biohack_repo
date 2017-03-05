# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi

# In demographic_models.py, we've defined a custom model for this problem
import demographic_models

import math

GreatGoodness = 20285641.5921
BestPopt = [1.881, 0.071, 1.845, 0.911, 0.355, 0.111]
iterrr = 0

data = dadi.Spectrum.from_file('YRI_CEU.fs')
ns = data.sample_sizes
pts_l = [40,50,60]
func = dadi.Demographics2D.split_mig
func = demographic_models.prior_onegrow_mig
func_ex = dadi.Numerics.make_extrap_log_func(func)

while (iterrr<len(BestPopt)):

	goodness = 0
	iterrr += 1
	popt = BestPopt
	popt[iterrr] += 0.05
	# Calculate the best-fit model AFS.
	model = func_ex(popt, ns, pts_l)
	# Likelihood of the data given the model AFS.
	ll_model = dadi.Inference.ll_multinom(model, data)
	# import math
	popt_nomig = array([ 1.897,  0.0388,  9.677,  0.395,  0.070])
	model_nomig = func_ex_nomig(popt_nomig, ns, pts_l)
	print(model_nomig)


	for i in xrange(len(model)):
		for j in xrange(len(model)):
			if math.isnan(data[i][j]):
				ddd = 0
			else:
				ddd = data[i][j]
			if math.isnan(model[i][j]):
				mmm = 0
			else:
				mmm = model[i][j]
			goodness += (ddd - mmm)*(ddd - mmm)
	
	if GreatGoodness < goodness:
		GreatGoodness = goodness
		BestPopt = popt
		iterrr = 0
		print GreatGoodness
		print BestPopt		
		continue

	goodness = 0
	popt = BestPopt
	popt[iterrr] -= 0.05
	# Calculate the best-fit model AFS.
	model = func_ex(popt, ns, pts_l)
	# Likelihood of the data given the model AFS.
	ll_model = dadi.Inference.ll_multinom(model, data)

	# The optimal value of theta given the model.

	# import math

	for i in xrange(len(model)):
		for j in xrange(len(model)):
			if math.isnan(data[i][j]):
				ddd = 0
			else:
				ddd = data[i][j]
			if math.isnan(model[i][j]):
				mmm = 0
			else:
				mmm = model[i][j]
			goodness += (ddd - mmm)*(ddd - mmm)
	
	if GreatGoodness < goodness:
		GreatGoodness = goodness
		BestPopt = popt
		iterrr = 0
		print GreatGoodness
		print BestPopt
		continue
	iterrr += 1


print(GreatGoodness)
print(BestPopt)

