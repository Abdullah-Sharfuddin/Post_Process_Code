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


# Read data
N = 128

x = []
y = []
z = []
rho = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rc, rr
    K = i - 1
    if (filesize("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-30.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-30.00-$K")
        append!(x,rr[:,2])
        append!(y,rr[:,3])
        append!(z,rr[:,4])
        append!(rho,rr[:,6])
    end
end

nn = Int(size(x,1))
count_w = 0
count_ice = 0

for i = 1:nn
    global count_w, count_ice
	if rho[i] == 1000.0
		count_w = count_w + 1
	else
		count_ice = count_ice + 1
	end
end

X_w = zeros(count_w)
X_ice = zeros(count_ice)
count_1 = 0
count_2 = 0
for i = 1:nn
    global X_w, X_ice, count_1, count_2
	if rho[i] == 1000.0
        count_1 = count_1 + 1
		# X_w[count_1] = z[i]
		X_w[count_1] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
	else
        count_2 = count_2 + 1
		# X_ice[count_2] = z[i]
		X_ice[count_2] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
	end
end


#=
# Find PDF
array_1 = X_w[1:count_w]
num_bins = 200
(prob_dens_1, bin_mid_1) = PDF_array(array_1,count_w,num_bins)

array_2 = X_ice[1:count_ice]
num_bins = 20
(prob_dens_2, bin_mid_2) = PDF_array(array_2,count_ice,num_bins)

# xticks=0:3:18,xlims=[0,18]
# yticks=-0.5:1:6.5,ylims=[-0.5,6.5]
plot(bin_mid_1[:],log.(10,prob_dens_1[:]),color = :bwr,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(bin_mid_2[:],log.(10,prob_dens_2[:]),color = :purple,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
=#