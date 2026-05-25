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

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,8])
end

n = 256
scalar = zeros(n,n,n)

scalar = nodal_value(var,scalar,Nx,Ny,Nz)


# # Compute liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# scalar .= C .* scalar

# Compute gradients
L = 0.512
h = L/n

gradient_x = zeros(n-2,n-2,n-2)
gradient_y = zeros(n-2,n-2,n-2)
gradient_z = zeros(n-2,n-2,n-2)
gradient = zeros(n-2,n-2,n-2)
for k = 2:n-1
    k1 = k - 1
    k2 = k + 1
    for j = 2:n-1
        j1 = j - 1
        j2 = j + 1
        for i = 2:n-1
            global gradient
            i1 = i - 1
            i2 = i + 1
            gradient_x[i-1,j-1,k-1] = (scalar[i2,j,k] - scalar[i1,j,k])/(2*h)
            gradient_y[i-1,j-1,k-1] = (scalar[i,j2,k] - scalar[i,j1,k])/(2*h)
            gradient_z[i-1,j-1,k-1] = (scalar[i,j,k2] - scalar[i,j,k1])/(2*h)
            gradient[i-1,j-1,k-1] = sqrt(gradient_x[i-1,j-1,k-1]^2 + gradient_y[i-1,j-1,k-1]^2 + gradient_z[i-1,j-1,k-1]^2)
        end
    end
end

# Compute the mean
mean_gradient = sum(gradient_z[:,:,:]) / ((n-2)^3)


# Compute PDF
array = gradient_z[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -14:2:8
y_labels = ["10^{-14}", "10^{-12}", "10^{-10}", "10^{-8}"]

plot(bin_mid[:] ,y,color = :summer,linestyle =:solid,linewidth=2,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,8])
end

n = 256
scalar = zeros(n,n,n)

scalar = nodal_value(var,scalar,Nx,Ny,Nz)


# # Compute liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# scalar .= C .* scalar

# Compute gradients
L = 0.512
h = L/n

gradient_x = zeros(n-2,n-2,n-2)
gradient_y = zeros(n-2,n-2,n-2)
gradient_z = zeros(n-2,n-2,n-2)
gradient = zeros(n-2,n-2,n-2)
for k = 2:n-1
    k1 = k - 1
    k2 = k + 1
    for j = 2:n-1
        j1 = j - 1
        j2 = j + 1
        for i = 2:n-1
            global gradient
            i1 = i - 1
            i2 = i + 1
            gradient_x[i-1,j-1,k-1] = (scalar[i2,j,k] - scalar[i1,j,k])/(2*h)
            gradient_y[i-1,j-1,k-1] = (scalar[i,j2,k] - scalar[i,j1,k])/(2*h)
            gradient_z[i-1,j-1,k-1] = (scalar[i,j,k2] - scalar[i,j,k1])/(2*h)
            gradient[i-1,j-1,k-1] = sqrt(gradient_x[i-1,j-1,k-1]^2 + gradient_y[i-1,j-1,k-1]^2 + gradient_z[i-1,j-1,k-1]^2)
        end
    end
end

# Compute the mean
mean_gradient = sum(gradient_z[:,:,:]) / ((n-2)^3)


# Compute PDF
array = gradient_z[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -14:2:8
y_labels = ["10^{-14}", "10^{-12}", "10^{-10}", "10^{-8}"]

plot!(bin_mid[:] ,y,color = :bwr,linestyle =:dash,linewidth=2,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

#=
# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P0_256/nodal_values/nodal_values-8.00-$I")
    append!(var,VAR[:,6])
end

n = 256
scalar = zeros(n,n,n)

scalar = nodal_value(var,scalar,Nx,Ny,Nz)


# # Compute liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# scalar .= C .* scalar

# Compute gradients
L = 0.512
h = L/n

gradient_x = zeros(n-2,n-2,n-2)
gradient_y = zeros(n-2,n-2,n-2)
gradient_z = zeros(n-2,n-2,n-2)
gradient = zeros(n-2,n-2,n-2)
for k = 2:n-1
    k1 = k - 1
    k2 = k + 1
    for j = 2:n-1
        j1 = j - 1
        j2 = j + 1
        for i = 2:n-1
            global gradient
            i1 = i - 1
            i2 = i + 1
            gradient_x[i-1,j-1,k-1] = (scalar[i2,j,k] - scalar[i1,j,k])/(2*h)
            gradient_y[i-1,j-1,k-1] = (scalar[i,j2,k] - scalar[i,j1,k])/(2*h)
            gradient_z[i-1,j-1,k-1] = (scalar[i,j,k2] - scalar[i,j,k1])/(2*h)
            gradient[i-1,j-1,k-1] = sqrt(gradient_x[i-1,j-1,k-1]^2 + gradient_y[i-1,j-1,k-1]^2 + gradient_z[i-1,j-1,k-1]^2)
        end
    end
end

# Compute the mean
mean_gradient = sum(gradient_z[:,:,:]) / ((n-2)^3)


# Compute PDF
array = gradient_z[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -6:1:0
y_labels = ["10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}", "10^{-2}", "10^{-1}", "10^{0}"]

plot!(bin_mid[:] ,y,color = :purple,linestyle =:dot,linewidth=2.5,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
=#