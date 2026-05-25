using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra

#=
# Mixed-phase parameters at 4500 ft altitude
L_w = 2.50*1e6
L_ice = 2.83*1e6
Rv = 461.5
R = 287.0
k = 0.0238
Cp = 1005.0
rho_a = 0.7674
rho_w = 1000.0
rho_ice = 917.0
mu_v = 0.000022
mu = 1.64e-5
T = 259.53
p = 57160
L = 0.512

# Calculate parameters
nu = mu/rho_a
p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

qvs_w = 621.97*p_sat_w/(p-p_sat_w)
qvs_ice = 621.97*p_sat_ice/(p-p_sat_ice)

qv_cloud = qvs_w * 1.01
qv_env = 1.75

qv_m = (qv_cloud+qv_env)/2

Se_w = qv_m/qvs_w - 1
Se_ice = qv_m/qvs_ice - 1

# Se_w = qv_cloud/qvs_w - 1
# Se_ice = qv_cloud/qvs_ice - 1
=#


# Altitude dependency
T0 = 288.15
p0 = 101325
# z1 = 2677.44
z1 = 3657.6
z2 = 7315.2
B = 0.0065
R = 8.314
Rd = 287
g = 9.8          #gravity
mu0 = 1.789e-5
S = 110.4        #Sutherland's constant
D0 = 2.2e-5
cp = 1005
k = 0.0238

T1 = T0 - B*z1
T2 = T0 - B*z2

p1 = p0*(1-(B*z1/T0))^(g/(Rd*B))
p2 = p0*(1-(B*z2/T0))^(g/(Rd*B))

rho_a_1 = p1/(Rd*T1)
rho_a_2 = p2/(Rd*T2)

mu1 = mu0*((T0+S)/(T1+S))*(T1/T0)^(3/2)
mu2 = mu0*((T0+S)/(T2+S))*(T2/T0)^(3/2)

D1 = D0*(p0/p1)*(T1/T0)^1.8
D2 = D0*(p0/p2)*(T2/T0)^1.8

alpa1 = k/(rho_a_1*cp)
alpa2 = k/(rho_a_2*cp)

Lw1 = 2501 - 2.36*(T1-273.15)
Lw2 = 2501 - 2.36*(T2-273.15)

Lice1 = 2834.1 - 0.29*(T1-273.15)
Lice2 = 2834.1 - 0.29*(T2-273.15)

# Calculate parameters
nu = mu1/rho_a_1
p_sat_w = 611.2*exp(17.67*(T1-273.15)/(T1-29.65))
p_sat_ice = 611.2*exp(22.46*(T1-273.15)/(T1-0.53))

qvs_w = 621.97*p_sat_w/(p1-p_sat_w)
qvs_ice = 621.97*p_sat_ice/(p1-p_sat_ice)

qv_cloud = qvs_w * 1.02
qv_env = 2.5

qv_m = (qv_cloud+qv_env)/2

Se_w = qv_m/qvs_w - 1
Se_ice = qv_m/qvs_ice - 1



#=
# Read data

temperature1 = readdlm("PR_DNS_Mixed/F3_CaseLH_128/temperature")
vapor1 = readdlm("PR_DNS_Mixed/F3_CaseLH_128/vapor")
supersat_i1 = readdlm("PR_DNS_Mixed/F3_CaseLH_128/supersat_ice")
mixed_rad1 = readdlm("PR_DNS_Mixed/F3_CaseLH_128/mixed-phase-radii",skipstart=1)
mixed_count1 = readdlm("PR_DNS_Mixed/F3_CaseLH_128/mixed-phase-count",skipstart=1)

temperature2 = readdlm("PR_DNS_Mixed/F3_CaseHH_128/temperature")
vapor2 = readdlm("PR_DNS_Mixed/F3_CaseHH_128/vapor")
supersat_i2 = readdlm("PR_DNS_Mixed/F3_CaseHH_128/supersat_ice")
mixed_rad2 = readdlm("PR_DNS_Mixed/F3_CaseHH_128/mixed-phase-radii",skipstart=1)
mixed_count2 = readdlm("PR_DNS_Mixed/F3_CaseHH_128/mixed-phase-count",skipstart=1)

temperature3 = readdlm("PR_DNS_Mixed/F3_CaseLL_128/temperature")
vapor3 = readdlm("PR_DNS_Mixed/F3_CaseLL_128/vapor")
supersat_i3 = readdlm("PR_DNS_Mixed/F3_CaseLL_128/supersat_ice")
mixed_rad3 = readdlm("PR_DNS_Mixed/F3_CaseLL_128/mixed-phase-radii",skipstart=1)
mixed_count3 = readdlm("PR_DNS_Mixed/F3_CaseLL_128/mixed-phase-count",skipstart=1)

temperature4 = readdlm("PR_DNS_Mixed/F3_CaseHL_128/temperature")
vapor4 = readdlm("PR_DNS_Mixed/F3_CaseHL_128/vapor")
supersat_i4 = readdlm("PR_DNS_Mixed/F3_CaseHL_128/supersat_ice")
mixed_rad4 = readdlm("PR_DNS_Mixed/F3_CaseHL_128/mixed-phase-radii",skipstart=1)
mixed_count4 = readdlm("PR_DNS_Mixed/F3_CaseHL_128/mixed-phase-count",skipstart=1)



# Calculate saturation vapor pressures
qvs_w = zeros(106)
qvs_ice = zeros(106)
total_mass = zeros(106)

for i=1:106
    global T,p,p_sat_ice,G_ice,sf
    T = temperature3[i,2]
    m_w = mixed_count1[i,5]
    m_ice = mixed_count1[i,6]
    qv = vapor1[i,2]
    p = 57160
    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))
    qvs_w[i] = 621.97*p_sat_w/(p-p_sat_w)
    qvs_ice[i] = 621.97*p_sat_ice/(p-p_sat_ice)
    total_mass[i] = m_w + m_ice + (qv*(L^3) *0.7674 *1e-3)
end

# plot(temperature1[1:106,1],qvs_w[1:106],color = :orange3,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(temperature1[1:106,1],qvs_ice[1:106],color = :bwr,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(temperature1[1:106,1],(qvs_w[1:106] ./ qvs_ice[1:106]),color = :black,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

plot(temperature1[1:106,1],total_mass[1:106],color = :nuuk,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
=#

#=
# Calculate growth factor
G_ice = zeros(106)
dr_dt = zeros(106)

for i=1:106
    global T,p_sat_ice,G_ice,sf
    T = temperature3[i,2]
    S = supersat_i3[i,2]
    r = mixed_rad3[i,6]
    sf = 1.0                  # Shape factor for sphere
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))
    G = ((L_ice*rho_ice)/(k*T))*((L_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*sf*p_sat_ice)
    G_ice[i] = 1/G
    dr_dt[i] = (G_ice[i]/r)*S*r*r
end

plot(temperature1[1:106,1],log.(10,dr_dt),color = :red4,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


G_ice = zeros(106)
dr_dt = zeros(106)

for i=1:106
    global T,p_sat_ice,G_ice,sf
    T = temperature7[i,2]
    S = supersat_i7[i,2]
    r = mixed_rad7[i,6]
    sf = 1.148              # Shape factor for hexagonal column
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))
    G = ((L_ice*rho_ice)/(k*T))*((L_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*sf*p_sat_ice)
    G_ice[i] = 1/G
    dr_dt[i] = (G_ice[i]/r)*S*r*r
end

plot!(temperature5[1:106,1],log.(10,dr_dt),color = :black,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


G_ice = zeros(106)
dr_dt = zeros(106)

for i=1:106
    global T,p_sat_ice,G_ice,sf
    T = temperature8[i,2]
    S = supersat_i8[i,2]
    r = mixed_rad8[i,6]
    sf = 1.082              # Shape factor for hexagonal plate
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))
    G = ((L_ice*rho_ice)/(k*T))*((L_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*sf*p_sat_ice)
    G_ice[i] = 1/G
    dr_dt[i] = (G_ice[i]/r)*S*r*r
end

plot!(temperature6[1:106,1],log.(10,dr_dt),color = :nuuk,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
=#


#=
P1 = readdlm("Post_Process_Scalar/validation/Chen_ri.txt",skipstart=1)
P2 = readdlm("Post_Process_Scalar/validation/mixed_ri.txt",skipstart=1)

plot(P1[:,1],P1[:,2],color = :red4,linestyle =:solid,linewidth=2,
labels="Chen,2024",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(P2[:,1],P2[:,2],color = :blue4,linestyle =:dash,linewidth=2,yticks=0:2:12,ylims=[0,12],
labels="Sharfuddin,2025",grid=false,thickness_scaling=1.5,legend=false)
# plot!(transition3[:,1],transition3[:,3],color = :green4,linestyle =:dot,linewidth=2.5,
# labels="L-L",grid=false,thickness_scaling=1.5,legend=false)
# plot!(transition4[:,1],transition4[:,3],color = :black,linestyle =:dashdot,linewidth=2,
# labels="H-L",grid=false,thickness_scaling=1.5)
=#

#=
# Validation
# IWC = 0.0086 @ 40 s

# yticks=8.6:0.2:10,ylims=[8.6,10]
A1 = readdlm("Post_Process_Scalar/Validation/Mixed/ri.txt",skipstart=1)
A2 = readdlm("Post_Process_Scalar/Validation/Mixed/n_ri.txt",skipstart=1)

A1[1,2] = 1
A2[1,2] = 1

plot(A1[:,1],A1[:,2],color = :red4,linestyle =:dash,linewidth=2,legend=false,
labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(A2[:,1],A2[:,2],color = :blue4,linestyle =:dot,linewidth=2.5,
labels="false",grid=false,thickness_scaling=1.5,yticks=0:2:13,ylims=[0,13])
=#





