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

#=
# Thermodynamic Forcing
# Read data
N = 128

# rd = []
# N1 = Int(N/4) + 1
# N2 = Int(N*3/4)
# for m = N1:N2
#     global k, rd
#     k = m - 1
#     rr = readdlm("PR_DNS_Scalar/Case_L_1/particle/particle_radius-0.00-$k")
#     append!(rd,rr[:,1])
# end

rc = []
r = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rc, r
    K = i - 1
    if (filesize("PR_DNS_Scalar/Case_L_2/particle/particle_radius-10.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Scalar/Case_L_2/particle/particle_radius-10.00-$K")
        append!(rc,rr[:,1])
        append!(r,rr[:,2])
    end
end

# num_drops = size(rd,1)
# rrc = zeros(num_drops)
# rra = zeros(num_drops)
# countc = 0
# counta = 0

# for i = 1:num_drops
#     global rrc, rra, countc, counta
    
# 	if r[i] > rc[i]
# 		countc = countc + 1
# 		rrc[countc] = r[i] 
# 	end
# 	if r[i] <= rc[i] && r[i] > rd[i]
# 		counta = counta + 1
# 		rra[counta] = r[i] 
# 	end
# 	if r[i] <= rd[i]
# 		counta = counta + 1
# 		rra[counta] = rd[i]
# 	end
# end


# Find PDF
# arrayc = rrc[1:countc]
# num_bins = 300
# (prob_dens_c, bin_mid_c) = PDF_array(arrayc,countc,num_bins)

# arraya = rra[1:counta]
# num_bins = 300
# (prob_dens_a, bin_mid_a) = PDF_array(arraya,counta,num_bins)

array = r[:]
count = size(array,1)
num_bins = 300
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)
V = 0.512^3
nc = (count/1e6)/V
# nc = 119
prob_dens .= prob_dens .*nc

# Find Mean
# meanc = sum(arrayc[:])/countc
# meana = sum(arraya[:])/counta

mean = sum(array[:])/count


# xticks=0:3:18,xlims=[0,18]
# yticks=-0.5:1:6.5,ylims=[-0.5,6.5]

plot(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :bwr,linestyle =:solid,linewidth=2,
labels="Cloud",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))

# plot!(bin_mid_a[:].*1e6,log.(10,prob_dens_a[:]),color = :purple,linestyle =:dash,linewidth=2,
# labels="Aerosol",grid=false,thickness_scaling=1.5,yticks=-1:2:7,ylims=[-1,7])
# plot!(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :orange3,linestyle =:dot,linewidth=2.5,
# labels="Critical",grid=false,thickness_scaling=1.5,xticks=0:2:12,xlims=[0,12])
# scatter!(bin_mid[1:10:end].*1e6,log.(10,prob_dens[1:10:end]),color = :orange,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# vline!([meanc*1e6],color = :summer,linestyle =:dot,linewidth=2.5,
# labels="Mean Cloud",grid=false,thickness_scaling=1.5,legend=false)
# vline!([meana*1e6],color = :black,linestyle =:dash,linewidth=2,
# labels="Mean Aerosol",grid=false,thickness_scaling=1.5)

x1 = [15, 15]; y1 = [3,8]
x2 = [1.26, 1.26]; y2 = [3,8]

plot!(x1,y1, color = :yellow4,linestyle =:dot,linewidth=2.5,
labels="H-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(x2,y2,color = :purple,linestyle =:dash,linewidth=2,
labels="H-H",grid=false,thickness_scaling=1.5)


r2 = []
N1 = 1
N2 = N
for i = N1:N2
    global K, r2
    K = i - 1
    if (filesize("PR_DNS_Scalar/Case_L_2/particle/particle_radius-25.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Scalar/Case_L_2/particle/particle_radius-25.00-$K")
        append!(r2,rr[:,2])
    end
end

array = r2[:]
count = size(array,1)
num_bins = 300
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)
V = 0.512^3
nc = (count/1e6)/V
# nc = 119
prob_dens .= prob_dens .*nc

# Find Mean
mean = sum(array[:])/count

plot!(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :summer,linestyle =:dashdot,linewidth=2,
labels="Cloud",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


r3 = []
N1 = 1
N2 = N
for i = N1:N2
    global K, r3
    K = i - 1
    if (filesize("PR_DNS_Scalar/Case_L_2/particle/particle_radius-35.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Scalar/Case_L_2/particle/particle_radius-35.00-$K")
        append!(r3,rr[:,2])
    end
end

array = r3[:]
count = size(array,1)
num_bins = 300
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)
V = 0.512^3
nc = (count/1e6)/V
# nc = 119
prob_dens .= prob_dens .*nc

# Find Mean
mean = sum(array[:])/count

plot!(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :black,linestyle =:dashdot,linewidth=2,
grid=false,legend=false)
scatter!(bin_mid[1:10:end].*1e6,log.(10,prob_dens[1:10:end]),color = :black,markershape = :circle,markersize = :2,
grid=false,legend=false,xticks=0:3:15,xlims=[-0.1,16])
# xticks=0:3:15,xlims=[-0.1,15.5]
# xticks=0:1:6,xlims=[0,6]
=#

#=
num_drops = size(rd,1)
array = rd[1:num_drops]
num_bins = 400
(prob_dens, bin_mid) = PDF_array(array,num_drops,num_bins)
mean = sum(array[:])/num_drops

plot(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :purple,linestyle =:dashdot,linewidth=2,xticks=0:1:5.5,xlims=[0,5.5],
labels="Aerosol",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
vline!([mean*1e6],color = :black,linestyle =:dash,linewidth=2,
labels="Mean Aerosol",grid=false,thickness_scaling=1.5)

vline([10],color = :summer,linestyle =:dot,linewidth=2.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"),
labels="Mean Cloud",grid=false,thickness_scaling=1.5,xticks=0:2:11,xlims=[0,11])
=#



#=
# Mixed-Phase Cloud
# Read data
N = 128

r = []
rho = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rho, r, rr
    K = i - 1
    if (filesize("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-30.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-30.00-$K")
        append!(r,rr[:,1])
		append!(rho,rr[:,6])
    end
end

num_drops = size(r,1)
rw = zeros(num_drops)
ri = zeros(num_drops)
countw = 0
counti = 0

for i = 1:num_drops
    global rw, ri, countw, counti
    
	if rho[i] == 1000.0
		countw = countw + 1
		rw[countw] = r[i] 
	else
		counti = counti + 1
		ri[counti] = r[i] 
	end
end

# Standard deviation
mean_w = sum(rw[:])/countw
sigma_w = 0
for i = 1:countw
	global sigma_w
	sigma_w = sigma_w + (rw[i]-mean_w)^2
end
sigma_w = sqrt(sigma_w/countw)

mean_i = sum(ri[:])/counti
sigma_i = 0
for i = 1:counti
	global sigma_i
	sigma_i = sigma_i + (ri[i]-mean_i)^2
end
sigma_i = sqrt(sigma_i/counti)


# Find PDF
arrayw = rw[1:countw]
num_bins = 400
(prob_dens_w1, bin_mid_w1) = PDF_array(arrayw,countw,num_bins)

arrayi = ri[1:counti]
num_bins = 40
# num_bins = 80
(prob_dens_i1, bin_mid_i1) = PDF_array(arrayi,counti,num_bins)

# Find Mean
meanw1 = sum(arrayw[:])/countw
meani1 = sum(arrayi[:])/counti

plot(bin_mid_w1[:].*1e6,prob_dens_w1[:].*1e-5,color = :orange3,linestyle =:dash,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(bin_mid_i1[:].*1e6,prob_dens_i1[:].*1e-5,color = :summer,linestyle =:solid,linewidth=2,
labels="H-H",xticks=0:2:16,xlims=[-0.2,16.2])
# xticks=0:2:16,xlims=[-0.2,16],yticks=-0.5:1:7,ylims=[-0.5,7]
# log.(10,

# x1 = [meanw1*1e6, meanw1*1e6]; y1 = [1,2]
# x2 = [meani1*1e6, meani1*1e6]; y2 = [1,2]

# x3 = [8, 8]; y3 = [1,6]
# x4 = [1, 1]; y4 = [1,6]

x1 = [0.000002261*1e6, 0.000002261*1e6]; y1 = [0,3]
x2 = [meani1*1e6, meani1*1e6]; y2 = [0,3]

x3 = [8, 8]; y3 = [0,14]
x4 = [1, 1]; y4 = [0,14]

plot!(x1,y1, color = :orange3,linestyle =:dash,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(x2,y2,color = :summer,linestyle =:solid,linewidth=2,
labels="H-H",grid=false,thickness_scaling=1.5)

plot!(x3,y3, color = :orange3,linestyle =:dashdot,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(x4,y4,color = :summer,linestyle =:dot,linewidth=2.5,
labels="H-H",grid=false,thickness_scaling=1.5)
=#


# Upscale Transfer of Energy and Thermal Boundary
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
		append!(Se,rr[:,4])
    end
end


# Find PDF
array = r[:]
count = size(array,1)
num_bins = 300
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

# Find Mean
mean = sum(array[:])/count

plot(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :summer,linestyle =:solid,linewidth=2,
grid=false,legend=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x = [mean*1e6, mean*1e6]; y = [3,6]

plot!(x,y, color = :summer,linestyle =:solid,linewidth=2,
grid=false,thickness_scaling=1.5,legend=false)



N = 512
r = []
Se = []

N1 = 1
N2 = N
for i = N1:N2
    global K, r, Se
    K = i - 1
    if (filesize("PR_DNS_Upscale/Case_L_P_256/particle/particle_radius-10.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Upscale/Case_L_P_256/particle/particle_radius-10.00-$K")
        append!(r,rr[:,2])
		append!(Se,rr[:,4])
    end
end


# Find PDF
array = r[:]
count = size(array,1)
num_bins = 300
(prob_dens, bin_mid) = PDF_array(array,count,num_bins)

# Find Mean
mean = sum(array[:])/count

plot!(bin_mid[:].*1e6,log.(10,prob_dens[:]),color = :bwr,linestyle =:dash,linewidth=2,
grid=false,legend=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

x = [mean*1e6, mean*1e6]; y = [3,6]

plot!(x,y, color = :bwr,linestyle =:dash,linewidth=2,
grid=false,thickness_scaling=1.5,legend=false)
