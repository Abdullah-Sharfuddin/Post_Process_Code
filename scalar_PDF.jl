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

    n = 288
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

#=
# Scalar Forcing
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
    VAR = readdlm("PR_DNS_Scalar/Case_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,5])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms
rms1 = rms

# Compute PDF
array = scalar[:]
num_bins = 200
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot(bin_mid[:],prob_dens[:],color = :red4,linestyle =:dot,linewidth=2.5,
labels="L-MGF",grid=false,thickness_scaling=1.5,legend=false)


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
    VAR = readdlm("PR_DNS_Scalar/Case_L_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,5])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms
rms2 = rms

# Compute PDF
array = scalar[:]
num_bins = 200
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot!(bin_mid[:],prob_dens[:],color = :blue4,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
# xticks=-6:2:4,xlims=[-6.1,4],yticks=-5:1:0,ylims=[-5,0]

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
    VAR = readdlm("PR_DNS_Scalar/Case_L_H_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,5])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end


# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms
rms3 = rms

# Compute PDF
array = scalar[:]
num_bins = 200
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot!(bin_mid[:],prob_dens[:],color = :green4,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


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
    VAR = readdlm("PR_DNS_Scalar/Case_H_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,5])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms
rms4 = rms

# Compute PDF
array = scalar[:]
num_bins = 200
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot!(bin_mid[:],prob_dens[:],color = :orange3,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)



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
    VAR = readdlm("PR_DNS_Scalar/Case_H_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,5])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms
rms5 = rms

# Compute PDF
array_5 = scalar[:]
num_bins = 200
count_5 = size(array_5,1)
(prob_dens_5, bin_mid_5) = PDF_array(array_5,count_5,num_bins)

plot!(bin_mid_5[:],prob_dens_5[:],color = :black,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
scatter!(bin_mid_5[1:6:end],prob_dens_5[1:6:end],color = :black,markershape = :square,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)
=#



# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
I = 0
for i = 1:N
    global I, var1, var2
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P2_288/nodal_values/nodal_values-6.50-$I")
    append!(var1,VAR[:,9])
    append!(var2,VAR[:,3])
end

n = 288
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var1,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:]) / (n*n*n)

# Compute PDF
array = scalar1[:]
num_bins = 400
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

plot(bin_mid[:] .*1e6,log.(10,prob_dens[:]),color = :summer,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x1 = [mean_scalar .*1e6, mean_scalar .*1e6]; y1 = [6,8]

plot!(x1,y1, color = :purple,linestyle =:dot,linewidth=2.5,
grid=false,thickness_scaling=1.5,legend=false)

