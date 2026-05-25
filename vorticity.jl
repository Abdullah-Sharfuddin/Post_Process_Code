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
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel1 = zeros(n,n,n)
vel2 = zeros(n,n,n)
vel3 = zeros(n,n,n)

vel_x = nodal_value(var1,vel1,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel2,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel3,Nx,Ny,Nz)

L = 0.512
h = L/n

# Compute velocity gradients
grad_vel_xx = zeros(n-1,n-1,n-1)
grad_vel_xy = zeros(n-1,n-1,n-1)
grad_vel_xz = zeros(n-1,n-1,n-1)
grad_vel_yx = zeros(n-1,n-1,n-1)
grad_vel_yy = zeros(n-1,n-1,n-1)
grad_vel_yz = zeros(n-1,n-1,n-1)
grad_vel_zx = zeros(n-1,n-1,n-1)
grad_vel_zy = zeros(n-1,n-1,n-1)
grad_vel_zz = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            grad_vel_xx[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xy[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xz[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_yx[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yy[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yz[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_zx[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zy[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zz[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
        end
    end
end

# Compute vorticity
vort = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            vort[i,j,k] = sqrt((grad_vel_zy[i,j,k] - grad_vel_yz[i,j,k])^2 + (grad_vel_xz[i,j,k] - grad_vel_zx[i,j,k])^2 + (grad_vel_yx[i,j,k] - grad_vel_xy[i,j,k])^2)
        end
    end
end

mean_vort= sum(vort[:,:,:]) / ((n-1)^3)

# Compute PDF
array = vort[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -6:1:-1
y_labels = ["10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}", "10^{-2}", "10^{-1}"]

plot(bin_mid[:],y,color = :summer,linestyle =:solid,linewidth=2,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x1 = [mean_vort, mean_vort]; y1 = [-6.5,-5]

plot!(x1,y1, color = :summer,linestyle =:solid,linewidth=2,
grid=false,thickness_scaling=1.5,legend=false)



# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel1 = zeros(n,n,n)
vel2 = zeros(n,n,n)
vel3 = zeros(n,n,n)

vel_x = nodal_value(var1,vel1,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel2,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel3,Nx,Ny,Nz)

L = 0.512
h = L/n

# Compute velocity gradients
grad_vel_xx = zeros(n-1,n-1,n-1)
grad_vel_xy = zeros(n-1,n-1,n-1)
grad_vel_xz = zeros(n-1,n-1,n-1)
grad_vel_yx = zeros(n-1,n-1,n-1)
grad_vel_yy = zeros(n-1,n-1,n-1)
grad_vel_yz = zeros(n-1,n-1,n-1)
grad_vel_zx = zeros(n-1,n-1,n-1)
grad_vel_zy = zeros(n-1,n-1,n-1)
grad_vel_zz = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            grad_vel_xx[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xy[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xz[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_yx[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yy[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yz[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_zx[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zy[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zz[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
        end
    end
end

# Compute vorticity
vort = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            vort[i,j,k] = sqrt((grad_vel_zy[i,j,k] - grad_vel_yz[i,j,k])^2 + (grad_vel_xz[i,j,k] - grad_vel_zx[i,j,k])^2 + (grad_vel_yx[i,j,k] - grad_vel_xy[i,j,k])^2)
        end
    end
end

mean_vort= sum(vort[:,:,:]) / ((n-1)^3)

# Compute PDF
array = vort[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -6:1:-1
y_labels = ["10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}", "10^{-2}", "10^{-1}"]

plot!(bin_mid[:],y,color = :bwr,linestyle =:dash,linewidth=2,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x1 = [mean_vort, mean_vort]; y1 = [-6.5,-5]

plot!(x1,y1, color = :bwr,linestyle =:dash,linewidth=2,
grid=false,thickness_scaling=1.5,legend=false)



# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P0_256/nodal_values/nodal_values-8.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel1 = zeros(n,n,n)
vel2 = zeros(n,n,n)
vel3 = zeros(n,n,n)

vel_x = nodal_value(var1,vel1,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel2,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel3,Nx,Ny,Nz)

L = 0.512
h = L/n


# Compute velocity gradients
grad_vel_xx = zeros(n-1,n-1,n-1)
grad_vel_xy = zeros(n-1,n-1,n-1)
grad_vel_xz = zeros(n-1,n-1,n-1)
grad_vel_yx = zeros(n-1,n-1,n-1)
grad_vel_yy = zeros(n-1,n-1,n-1)
grad_vel_yz = zeros(n-1,n-1,n-1)
grad_vel_zx = zeros(n-1,n-1,n-1)
grad_vel_zy = zeros(n-1,n-1,n-1)
grad_vel_zz = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            grad_vel_xx[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xy[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_xz[i,j,k] = (vel_x[i+1,j,k] - vel_x[i,j,k])/h
            grad_vel_yx[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yy[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_yz[i,j,k] = (vel_y[i,j+1,k] - vel_y[i,j,k])/h
            grad_vel_zx[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zy[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
            grad_vel_zz[i,j,k] = (vel_z[i,j,k+1] - vel_z[i,j,k])/h
        end
    end
end

# Compute vorticity
vort = zeros(n-1,n-1,n-1)

for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            vort[i,j,k] = sqrt((grad_vel_zy[i,j,k] - grad_vel_yz[i,j,k])^2 + (grad_vel_xz[i,j,k] - grad_vel_zx[i,j,k])^2 + (grad_vel_yx[i,j,k] - grad_vel_xy[i,j,k])^2)
        end
    end
end

mean_vort= sum(vort[:,:,:]) / ((n-1)^3)

# Compute PDF
array = vort[:]
num_bins = 300
count = size(array,1)
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

y = log.(10,prob_dens[:])
y_ticks_pre = -6:1:-1
y_labels = ["10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}", "10^{-2}", "10^{-1}"]

plot!(bin_mid[:],y,color = :purple,linestyle =:dot,linewidth=2.5,yticks = (y_ticks_pre, y_labels),
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x1 = [mean_vort, mean_vort]; y1 = [-6.5,-5]

plot!(x1,y1, color = :purple,linestyle =:dot,linewidth=2.5,
grid=false,thickness_scaling=1.5,legend=false)
