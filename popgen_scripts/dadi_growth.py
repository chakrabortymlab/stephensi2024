# growth: Exponential growth beginning some time ago.
# script for running dadi to fit a growth demographic model 
# to the SFS of synonymous SNPs from the Trivandrum population.

#Author:Alex Samano 2024

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D

#import sfs
filename= "/Users/Alex/Desktop/ChakLab/asteph_popgen/dadi_demography/tvm_synon.fs"
fs=dadi.Spectrum.from_file(filename) 

#mutation rate for An. stephensi from Rashid 2022
mu=1.36e-9
#since we are using the synonymous SFS, L=total length of coding sites/3
L=205167228

#sample size from sfs
ns = fs.sample_sizes
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] #define grid points

#define model
func = Demographics1D.growth
param_names= ("nu","T")
upper_bound = [10, 1] #set boundaries
lower_bound = [1e-4, 1e-5]
p0 = [0.01,0.001]
maxiter=100

# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
#perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

# optimize:
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=maxiter)

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data=dadi.Inference.ll_multinom(fs, fs)
# calculate best fit synonymous theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

# Model specific scaling of parameters 
Nanc=theta / (4*mu*L)
nu_scaled_dip=popt[0]*Nanc
T_scaled_gen=popt[1]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
scaled_popt=(Nanc,nu_scaled_dip,T_scaled_gen)

#print and plot residuals
print("ll",ll_model,popt,"scaled",scaled_popt)
dadi.Plotting.plot_1d_comp_multinom(model, fs)
