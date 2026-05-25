
using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra


# Given values at 4 s
TKE_L = 0.0025
TKE_H = 0.0094
eps_L = 0.00063
eps_H = 0.00723
eta_L = 0.00157
eta_H = 0.00091
urms_L = 0.0797
urms_H = 0.1386


# Stokes number
L = 0.512
rho_l = 1000
rho_a = 1
nu = 1.5e-5
g = 9.8
r = 10*1e-6       # mean value
urms_0 = 0.16        # initial rms velocity

tau_p = (2*rho_l*(r^2))/(9*rho_a*nu)
tau_f = L/urms_0

St = tau_p/tau_f

# For grid resolution
N = 256
kmax = (N*pi)/L
kmax_eta_L = kmax*eta_L
kmax_eta_H = kmax*eta_H

# Taylor microscale
lambda_L = (10*nu*TKE_L)/eps_L
lambda_L = sqrt(lambda_L)
lambda_H = (10*nu*TKE_H)/eps_H
lambda_H = sqrt(lambda_H)

# Integral scale
l_L = (TKE_L^(3/2))/eps_L
l_H = (TKE_H^(3/2))/eps_H

# Taylor Reynolds number
Re_lambda_L = (urms_L*lambda_L)/nu
Re_lambda_H = (urms_H*lambda_H)/nu

# Integral Reynolds number
Re_l_L = (urms_L*l_L)/nu
Re_l_H = (urms_H*l_H)/nu

# Forcing Reynolds number
k0 = (2*pi)/L
Re_f_L = eps_L^(1/3)*(k0^(-4/3))/nu
Re_f_H = eps_H^(1/3)*(k0^(-4/3))/nu




#=
# Check
# TKE = 15.586
# eps = 16.058
# eta = 0.0314
# urms = 5.58

TKE = 58.89
eps = 132.3
eta = 0.0175
urms = 6.266

# For grid resolution
N = 128
L = 2*pi
k0 = (2*pi)/L
kmax = (N*pi)/L
kmax_eta = kmax*eta

# Taylor microscale
nu = 2.5e-2
lambda = (10*nu*TKE)/eps
lambda = sqrt(lambda)

# Integral scale
l = (TKE^(3/2))/eps

# Taylor Reynolds number
Re_lambda = (urms*lambda)/nu

# Integral Reynolds number
Re_l = (urms*l)/nu

# Turbulent length scale
u_prime = (2*TKE/3)^(1/2)
L_eps = (u_prime^3)/eps

# Eddy turnover time
t_E = l/u_prime
=#

#=
# Normalized variance
beta = 1.0
temperature = readdlm("PR_DNS_Scalar/F3_CaseVAL_128/temperature", skipstart=2)
var = temperature[11,3]
var_normal = var/((beta*L_eps)^2)
=#

#=
L = 0.512
nu = 1.5e-5
urms_V = 0.1817
urms_H = 0.1239
urms_M = 0.0890
urms_L = 0.0605

Re_V = (urms_V*L)/nu
Re_H = (urms_H*L)/nu
Re_M = (urms_M*L)/nu
Re_L = (urms_L*L)/nu
=#

#=
# Check
nu = 1.5e-5
TKE_H = 0.0165
TKE_L = 0.00183
eps_H = 0.01172
eps_L = 0.00019
eta_H = 0.00073
eta_L = 0.002

# Large-eddy lengthscale
L_H = TKE_H^(3/2)/eps_H
L_L = TKE_L^(3/2)/eps_L

# Ratio
r_H = eta_H/L_H
r_L = eta_L/L_L

# Reynolds number
Re_H = TKE_H^2/(eps_H*nu)
Re_L = TKE_L^2/(eps_L*nu)
=#

