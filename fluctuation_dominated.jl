using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

# Function to compute PDF
function PDF_array(array,num,num_bins)

	var_min = 10^10
	var_max = -10^10

	array_size = num

	for i = 1:array_size
		if array[i] < var_min
			var_min = array[i]
		end
		if array[i] > var_max
			var_max = array[i]
		end
	end

	bin_size = (var_max - var_min)/(num_bins-1)
	
	prob_dens = zeros(num_bins,1)

	for i = 1:array_size
		for j = 1:num_bins
			if (array[i] >= var_min+(j-0.5)*bin_size) && (array[i] < var_min+(j+0.5)*bin_size)
				prob_dens[j] = prob_dens[j] + 1.0
				break
			end
		end
	end

	total_num = 0
	# Normalize Probability density function
	for j = 1:num_bins
		total_num = total_num + prob_dens[j]
	end

	for j = 1:num_bins
		prob_dens[j] = prob_dens[j]/(bin_size*total_num)
	end

	bin_mid = zeros(num_bins,1)

	for i = 1:num_bins
        bin_mid[i] = var_min+(0.5+i)*bin_size
    end

	return prob_dens, bin_mid
end

# Read nodal values
function nodal_value(var,scalar,Nx,Ny,Nz)

    n = 192
    nx = Int(n/Nx)
    ny = Int(n/Ny)
    nz = Int(n/Nz)

    count = 0

    for kk = 1:Nz
        k1 = 1 + nz*(kk-1)
        k2 = nz + nz*(kk-1)
        for jj = 1:Ny
            j1 = 1 + ny*(jj-1)
            j2 = ny + ny*(jj-1)
            for ii = 1:Nx
                i1 = 1 + nx*(ii-1)
                i2 = nx + nx*(ii-1)
                for k = k1:k2
                    for j = j1:j2
                        for i = i1:i2
                            count = count + 1
                            scalar[i,j,k] = var[count]
                        end
                    end
                end
            end
        end
    end

    return scalar, count
end

## PDF of environmental supersaturation
# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_L_2/nodal_values/nodal_values-5.00-$I")
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuation
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#         end
#     end
# end

# Compute PDF
array = scalar[:]
num_bins = 200
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot(bin_mid[:],prob_dens[:],color = :summer,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,legend=false)

x1 = [mean_scalar, mean_scalar]; y1 = [-4,-3]

plot!(x1,y1, color = :summer,linestyle =:solid,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,legend=false)

## Critical equilibrium supersaturation
sigma_l = 0.072
M_l = 0.018
R = 8.314
rho_l = 1000

temperature1 = readdlm("PR_DNS_Scalar/Case_H_L_2/temperature")
T = temperature1[11,2]

A = (2*sigma_l*M_l)/(R*T*rho_l)

rd = 0.1*1e-6
kappa = 0.61
Sc = (2/sqrt(kappa))*(A/(3*rd))^(3/2)

vline!([Sc],color = :purple,linestyle =:dot,linewidth=2.5,
labels="H-H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


#=
## PDF of equilibrium supersaturation

Sk = []
N1 = 1
N2 = N
for i = N1:N2
    global K, Sk
    K = i - 1
    if (filesize("PR_DNS_Scalar/Case_H_L_2/particle/particle_radius-25.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Scalar/Case_H_L_2/particle/particle_radius-25.00-$K")
        append!(Sk,rr[:,3])
    end
end

# Compute mean
mean = sum(Sk[:]) / (16e6)

array = Sk[:]
count = size(array,1)
num_bins = 200
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot!(bin_mid[:],log.(10,prob_dens[:]),color = :purple,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
=#


