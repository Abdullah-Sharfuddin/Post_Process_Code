
using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra

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

    n = 128
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

# instantaneous environmental supersaturation
N = 128
Nx = 8
Ny = 4
Nz = 4

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseHL1_2_128/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,6])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# instantaneous particle equilibrium supersaturation
Sk = []
N1 = 1
N2 = N
for i = N1:N2
    global I, rc, r
    I = i - 1
    if (filesize("PR_DNS_Scalar/F3_CaseHL1_2_128/particle/particle_radius-10.00-$I") == 0)
        continue
    else
        ss = readdlm("PR_DNS_Scalar/F3_CaseHL1_2_128/particle/particle_radius-10.00-$I")
        append!(Sk,ss[:,3])
    end
end

# Find mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_Sk = sum(Sk[:]) / size(Sk,1)

# Find PDF
count1 = n^3
array1 = scalar[:]
num_bins = 400
(prob_dens1, bin_mid1) = PDF_array(array1,count1,num_bins)

count2 = size(Sk,1)
array2 = Sk[:]
num_bins = 400
(prob_dens2, bin_mid2) = PDF_array(array2,count2,num_bins)


plot(bin_mid1[:],log.(10,prob_dens1[:]),color = :bwr,linestyle =:dashdot,linewidth=2,
labels="Se",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(bin_mid2[:],log.(10,prob_dens2[:]),color = :summer,linestyle =:dot,linewidth=2.5,
labels="Sk",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# xticks=-0.03:0.005:-0.01,xlims=[-0.03,-0.01]
vline!([mean_scalar],color = :black,linestyle =:dash,linewidth=2,
labels="Mean Se",grid=false,thickness_scaling=1.5,legend=false)
vline!([mean_Sk],color = :brown,linestyle =:solid,linewidth=2,
labels="Mean Sk",grid=false,thickness_scaling=1.5)


#=
plot(bin_mid1[:],log.(10,prob_dens1[:]),color = :summer,linestyle =:solid,linewidth=2,
labels="Se",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(bin_mid2[:],log.(10,prob_dens2[:]),color = :brown,linestyle =:dash,linewidth=2,
labels="Sk",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# xticks=-0.03:0.005:-0.01,xlims=[-0.03,-0.01]
vline!([mean_scalar],color = :bwr,linestyle =:dot,linewidth=2.5,
labels="Mean Se",grid=false,thickness_scaling=1.5,legend=false)
vline!([mean_Sk],color = :black,linestyle =:dashdot,linewidth=2,
labels="Mean Sk",grid=false,thickness_scaling=1.5)
=#