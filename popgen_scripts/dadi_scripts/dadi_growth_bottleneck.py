"""period of growth followed by an an instantaneous size change some time ago.
script for running dadi to fit a growth+bottleneck demographic model 
to the SFS of synonymous SNPs from the Lakshadweep population. This is the best fitting 1D model for LKD. Parameters were used for neutral simulations to estimate a significance threshold for G12.
Author: Alex Samano 2024"""


import dadi
import numpy
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D


#import sfs
filename= "lkd60_synon_proj36.fs"
fs=dadi.Spectrum.from_file(filename) 

#mutation rate for An. stephensi from Rashid 2022
mu=1.36e-9
#since we are using the synonymous SFS, L=total length of coding sites/3
L=17973548

ns = fs.sample_sizes
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]

def growth_bottleneck(params, ns, pts):
    """
    Exponential growth beginning some time ago.
    ns = (n1,)
    nuG: Ratio of growth to ancient population size
    nuB: Ratio of bottleneck pop size to ancient
    nuP: Ratio of present to ancient size
    Tg: Time in the past at which growth began (in units of 2*Na generations) 
    Tb: Length of bottleneck
    Tp: Time since bottleneck recovery

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuG, nuB , nuP, Tg, Tb, Tp= params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    #growth
    nu_func = lambda t: numpy.exp(numpy.log(nuG) * t/Tg)
    phi = Integration.one_pop(phi, xx, Tg, nu_func)
    #bottleneck
    phi = Integration.one_pop(phi, xx, Tb, nuB)
    #recovery
    phi = Integration.one_pop(phi, xx, Tp, nuP)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

func=growth_bottleneck
#bounds
upper_bound = [10,10,10,1,1, 1]
lower_bound = [1e-4,1e-4,1e-4,1e-5,1e-5,1e-5]
p0 = [1,0.1 ,0.1,0.1,0.07,0.01]
maxiter=100
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound,upper_bound=upper_bound, maxiter=maxiter)

model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data=dadi.Inference.ll_multinom(fs, fs)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)


# Model specific scaling of parameters (will depend on mu and L that you supply)
Nanc=theta / (4*mu*L)
nuG_scaled=popt[0]*Nanc
nuB_scaled=popt[1]*Nanc
nuP_scaled=popt[2]*Nanc
Tg_scaled=popt[3]*2*Nanc
Tb_scaled=popt[4]*2*Nanc
Tp_scaled=popt[5]*2*Nanc

scaled_popt=(Nanc,nuG_scaled,nuB_scaled,nuP_scaled,Tg_scaled,Tb_scaled,Tp_scaled)

#print scaled parameters and plot residuals
print("ll",ll_model,popt,"scaled",scaled_popt)
#dadi.Plotting.plot_1d_comp_multinom(model, fs)
