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

function nodal_value(var,scalar,Nx,Ny,Nz)

    n = 256
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

    return scalar
end


# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
var4 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,4])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)
scalar3 = zeros(n,n,n)
scalar4 = zeros(n,n,n)
scalar1 = nodal_value(var1,scalar1,Nx,Ny,Nz)
scalar2 = nodal_value(var2,scalar2,Nx,Ny,Nz)
scalar3 = nodal_value(var3,scalar3,Nx,Ny,Nz)
scalar4 = nodal_value(var4,scalar4,Nx,Ny,Nz)


# Compute mean field
mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)
mean_scalar3 = sum(scalar3[:,:,:]) / (n*n*n)
mean_scalar4 = sum(scalar4[:,:,:]) / (n*n*n)

T = mean_scalar4

# # Compute fluctuating field
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             global scalar1,scalar2,scalar3,scalar4
#             scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
#             scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
#             scalar3[i,j,k] = scalar3[i,j,k] - mean_scalar3
#             scalar4[i,j,k] = scalar4[i,j,k] - mean_scalar4
#         end
#     end
# end

# Compute TKE dissipation rate
L = 0.512
h = L/n

mu = 1.5e-5
rho_a = 1.0
nu = mu/rho_a

TKE_disp_rate = zeros(n-2,n-2,n-2)
comp_1 = zeros(n-2,n-2,n-2)
comp_2 = zeros(n-2,n-2,n-2)
comp_3 = zeros(n-2,n-2,n-2)
comp_4 = zeros(n-2,n-2,n-2)
comp_5 = zeros(n-2,n-2,n-2)
comp_6 = zeros(n-2,n-2,n-2)

for k = 2:n-1
    k1 = k - 1
    k2 = k + 1
    for j = 2:n-1
        j1 = j - 1
        j2 = j + 1
        for i = 2:n-1
            global TKE_disp_rate,comp_1,comp_2,comp_3,comp_4,comp_5,comp_6
            i1 = i - 1
            i2 = i + 1
            comp_1[i-1,j-1,k-1] = 2*((scalar1[i2,j,k] - scalar1[i1,j,k])/(2*h))^2
            comp_2[i-1,j-1,k-1] = 2*((scalar2[i,j2,k] - scalar2[i,j1,k])/(2*h))^2
            comp_3[i-1,j-1,k-1] = 2*((scalar3[i,j,k2] - scalar3[i,j,k1])/(2*h))^2
            comp_4[i-1,j-1,k-1] = ((scalar1[i,j2,k] - scalar1[i,j1,k])/(2*h) + (scalar2[i2,j,k] - scalar2[i1,j,k])/(2*h))^2
            comp_5[i-1,j-1,k-1] = ((scalar1[i,j,k2] - scalar1[i,j,k1])/(2*h) + (scalar3[i2,j,k] - scalar2[i1,j,k])/(2*h))^2
            comp_6[i-1,j-1,k-1] = ((scalar2[i,j,k2] - scalar2[i,j,k1])/(2*h) + (scalar3[i,j2,k] - scalar3[i,j1,k])/(2*h))^2
            TKE_disp_rate[i-1,j-1,k-1] = nu*(comp_1[i-1,j-1,k-1]+comp_2[i-1,j-1,k-1]+comp_3[i-1,j-1,k-1]+comp_4[i-1,j-1,k-1]+comp_5[i-1,j-1,k-1]+comp_6[i-1,j-1,k-1])
        end
    end
end

# Compute the timescales for large eddy turnover, Taylor microscale, and Kolmogorov scale
timescale_1 = zeros(n-2,n-2,n-2)
timescale_2 = zeros(n-2,n-2,n-2)
timescale_3 = zeros(n-2,n-2,n-2)

for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global timescale_1, timescale_2, timescale_3
            timescale_1[i-1,j-1,k-1] = (scalar1[i-1,j-1,k-1])^2 + (scalar2[i-1,j-1,k-1])^2 + (scalar3[i-1,j-1,k-1])^2
            timescale_1[i-1,j-1,k-1] = L/sqrt(timescale_1[i-1,j-1,k-1])
            timescale_2[i-1,j-1,k-1] = (15*nu)/TKE_disp_rate[i-1,j-1,k-1]
            timescale_2[i-1,j-1,k-1] = sqrt(timescale_2[i-1,j-1,k-1])
            timescale_3[i-1,j-1,k-1] = (nu/TKE_disp_rate[i-1,j-1,k-1])^0.5
        end
    end
end


# Compute the timescales for microscale phase change, particle response and Brownian diffusion
N = 1024
r = []
Se = []

N1 = 1
N2 = N
for i = N1:N2
    global K, r, Se
    K = i - 1
    if (filesize("PR_DNS_Upscale/Case_H_P_256/particle/particle_radius-10.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Upscale/Case_H_P_256/particle/particle_radius-10.00-$K")
        append!(r,rr[:,2])
        append!(Se,rr[:,2])
    end
end

nn = size(r,1)
timescale_4 = zeros(nn)
timescale_5 = zeros(nn)
timescale_6 = zeros(nn)

Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 2.2e-5
kB = 1.38e-23
nu = 1.5e-5
rho_a = 1.0

p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))
G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
G = 1/G

for i=1:nn
    global timescale_4, timescale_5, timescale_6
    timescale_4[i] = (r[i]^2)/abs(G*Se[i])
    timescale_5[i] = r[i]^3
    timescale_5[i] = timescale_5[i]*(6*pi*mu)/(kB*T)
    timescale_6[i] = (2*rho*r[i]^2)/(9*rho_a*nu)
end


# Compute PDF
array1 = timescale_1[:]
num_bins = 400
count1 = size(array1,1)
(prob_dens1, bin_mid1) = PDF_array(array1,count1,num_bins)

array2 = timescale_2[:]
num_bins = 400
count2 = size(array2,1)
(prob_dens2, bin_mid2) = PDF_array(array2,count2,num_bins)

array3 = timescale_3[:]
num_bins = 400
count3 = size(array3,1)
(prob_dens3, bin_mid3) = PDF_array(array3,count3,num_bins)

array4 = timescale_4[:]
num_bins = 400
count4 = size(array4,1)
(prob_dens4, bin_mid4) = PDF_array(array4,count4,num_bins)

array5 = timescale_5[:]
num_bins = 400
count5 = size(array5,1)
(prob_dens5, bin_mid5) = PDF_array(array5,count5,num_bins)

array6 = timescale_6[:]
num_bins = 400
count6 = size(array6,1)
(prob_dens6, bin_mid6) = PDF_array(array6,count6,num_bins)


x_ticks_pre = -4:2:6
x_labels = ["10^{-4}", "10^{-2}", "10^{0}", "10^{2}", "10^{4}", "10^{6}"]

y_ticks_pre = -10:4:6
y_labels = ["10^{-10}", "10^{-6}", "10^{-2}", "10^{2}", "10^{6}"]

x = log.(10,bin_mid1[:]); y = log.(10,prob_dens1[:])
plot(x,y,color = :purple,linestyle =:dot,linewidth=2.5,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x = log.(10,bin_mid2[:]); y = log.(10,prob_dens2[:])
# x = [array2[1], array2[1]]; y = [2,6]
plot!(x,y,color = :blue4,linestyle =:solid,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5)

x = log.(10,bin_mid3[:]); y = log.(10,prob_dens3[:])
# x = [array3[1], array3[1], array3[1]]; y = [2,4,6]
plot!(x,y,color = :orange3,linestyle =:dash,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
scatter!(x[1:40:end],y[1:40:end],color = :orange3,markershape = :circle,markersize = :2,
labels=false,grid=false,thickness_scaling=1.5,xlims=(-4,6.2))

x = log.(10,bin_mid4[:]); y = log.(10,prob_dens4[:])
plot!(x,y,color = :green4,linestyle =:dash,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5)

x = log.(10,bin_mid5[:]); y = log.(10,prob_dens5[:])
plot!(x,y,color = :red4,linestyle =:solid,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5)
scatter!(x[1:40:end],y[1:40:end],color = :red4,markershape = :square,markersize = :2,
labels=false,grid=false,thickness_scaling=1.5)

x = log.(10,bin_mid6[:]); y = log.(10,prob_dens6[:])
plot!(x,y,color = :black,linestyle =:dashdot,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5)


