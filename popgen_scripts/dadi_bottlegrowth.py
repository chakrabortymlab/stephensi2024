# script for running dadi to fit a bottlegrowth demographic model 
# to the SFS of synonymous SNPs from the Trivandrum population.
# Was run 100 times, parameters from run with highest likelihood
# were used for ABC simulations and estimating selection coefficients.


import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D

filename= "/asteph_popgen/dadi_demography/tvm_synon.fs"
fs=dadi.Spectrum.from_file(filename) 

#mutation rate for An. stephensi from Rashid 2022
mu=1.36e-9
#since we are using the synonymous SFS, L=total length of coding sites/3
L=17973548


ns = fs.sample_sizes #get sample size from data
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] #define grid points
func = Demographics1D.bottlegrowth 
param_names= ("nuB","nuF","T")  
upper_bound = [10, 10, 0.1] #define boundaries
lower_bound = [1e-4, 1e-4, 1e-5]
p0 = [0.01,0.1,0.005] # initial parameters

maxiter=100
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
#run optimization
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=maxiter)

model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data=dadi.Inference.ll_multinom(fs, fs)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

#scale params
Nanc=theta / (4*mu*L)
nuB_scaled_dip=popt[0]*Nanc
nuF_scaled_dip=popt[1]*Nanc
T_scaled_gen=popt[2]*2*Nanc  
scaled_popt=(Nanc,nuB_scaled_dip,nuF_scaled_dip,T_scaled_gen)
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
print("ll",ll_model,popt,"scaled",scaled_popt)

#create plot
#dadi.Plotting.plot_1d_comp_multinom(model, fs)


