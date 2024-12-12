# script for running dadi to fit a three-epoch demographic model 
# to the SFS of synonymous SNPs from the Trivandrum population.

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D

filename= "/asteph_popgen/dadi_demography/tvm_synon.fs"
fs=dadi.Spectrum.from_file(filename) 

#mutation rate for An. stephensi from Rashid 2022
mu=1.36e-9
#since we are using the synonymous SFS, L=total length of coding sites/3
L=205167228

ns = fs.sample_sizes
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] #define grid points
func = Demographics1D.three_epoch
param_names= ("nuB","nuF","TB","TF")
upper_bound = [10, 10, 0.1, 0.1] #set boundaries
lower_bound = [1e-4, 1e-4, 1e-5, 1e-5]
p0 = [0.01,0.1,0.005,0.001] # initial parameters

maxiter=100
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)
#perturb params
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
#run optimization
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=maxiter)

#calculate synonymous theta
model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data=dadi.Inference.ll_multinom(fs, fs)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

Nanc=theta / (4*mu*L)
nuB_scaled_dip=popt[0]*Nanc
nuF_scaled_dip=popt[1]*Nanc
TB_scaled_gen=popt[2]*2*Nanc  
TF_scaled_gen=popt[3]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuF_scaled_dip","TB_scaled_gen","TF_scaled_gen")
scaled_popt=(Nanc,nuB_scaled_dip,nuF_scaled_dip,TB_scaled_gen,TF_scaled_gen)
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
print("ll",ll_model,popt,"scaled",scaled_popt)
dadi.Plotting.plot_1d_comp_multinom(model, fs)




java -jar msms.jar msms 36 1 -N 1757732.40914425 -t 382.48257222979 -r 28123.7185463081 -s 753 -I 1 36 -en 0.0004881430665 1 0.00130724133 -eg 0.0004881430665 1 14941.0606628007 -Sp 0.5 -SI tbs 1 5.6891481023943e-07 -SFC -SAA tbs -SAa tbs -oTrace



msstats=/usr/local/bin/libsequence/examples/msstats