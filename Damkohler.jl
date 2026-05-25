using DelimitedFiles
using Plots
using LaTeXStrings
using LinearAlgebra
using CurveFit



#=
A1 = readdlm("Post_Process_Scalar/Data/Da_V2.txt",skipstart=1)
A2 = readdlm("Post_Process_Scalar/Data/Da_L2.txt",skipstart=1)

B1 = readdlm("Post_Process_Scalar/Data/urms_V1.txt",skipstart=1)
B2 = readdlm("Post_Process_Scalar/Data/urms_L1.txt",skipstart=1)

# n = 40
n = 38
Da1 = 10 .^A1[1:n,2]
Da2 = 10 .^A2[1:n,2]

L = 0.512
tt1 = L ./B1[1:n,2]
tt2 = L ./B2[1:n,2]

tc1 = tt1 ./ Da1
tc2 = tt2 ./ Da2

D1 = tt1 ./ tc1
D2 = tt2 ./ tc2

# plot(A1[1:n,1],log.(10,tt1),color = :red4,linestyle =:dot,linewidth=2.5,legendfont=font(5,"Helvetica Bold"),
# labels="τ_t, V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(A2[1:n,1],log.(10,tt2),color = :blue4,linestyle =:dash,linewidth=2,
# labels="τ_t, L-L",grid=false,thickness_scaling=1.5,xticks=0:5:20,xlims=[0,20])
plot(A1[1:n,1],log.(10,tc1),color = :green4,linestyle =:solid,linewidth=2,
labels="τ_c, V-L",grid=false,thickness_scaling=1.5,xticks=0:5:20,xlims=[0,20])
plot!(A2[1:n,1],log.(10,tc2),color = :black,linestyle =:dashdot,linewidth=2,
labels="τ_c, L-L",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"),tickfont=font(9,"Helvetica Bold"))

# plot(A1[:,1],D1,color = :blue4,linestyle =:dash,linewidth=2,
# labels="V",grid=false,thickness_scaling=1.5)
=#



#=
u1 = readdlm("Post_Process_Scalar/Data/urms_V1.txt",skipstart=1)
u2 = readdlm("Post_Process_Scalar/Data/urms_H1.txt",skipstart=1)
u3 = readdlm("Post_Process_Scalar/Data/urms_H1_prime.txt",skipstart=1)
u4 = readdlm("Post_Process_Scalar/Data/urms_M1_prime.txt",skipstart=1)
u5 = readdlm("Post_Process_Scalar/Data/urms_M1.txt",skipstart=1)
u6 = readdlm("Post_Process_Scalar/Data/urms_L1.txt",skipstart=1)


D1 = readdlm("Post_Process_Scalar/Data/Da_V2.txt",skipstart=1)
D2 = readdlm("Post_Process_Scalar/Data/Da_H2.txt",skipstart=1)
D3 = readdlm("Post_Process_Scalar/Data/Da_H2_prime.txt",skipstart=1)
D4 = readdlm("Post_Process_Scalar/Data/Da_M2_prime.txt",skipstart=1)
D5 = readdlm("Post_Process_Scalar/Data/Da_M2.txt",skipstart=1)
D6 = readdlm("Post_Process_Scalar/Data/Da_L2.txt",skipstart=1)

L = 0.512
tt1 = L ./u1[1:37,2]
tt2 = L ./u2[1:37,2]
tt3 = L ./u3[1:37,2]
tt4 = L ./u4[1:37,2]
tt5 = L ./u5[1:37,2]
tt6 = L ./u6[1:37,2]

tc1 = tt1 ./ (10 .^ D1[1:37,2])
tc2 = tt2 ./ (10 .^ D2[1:37,2])
tc3 = tt3 ./ (10 .^ D3[1:37,2])
tc4 = tt4 ./ (10 .^ D4[1:37,2])
tc5 = tt5 ./ (10 .^ D5[1:37,2])
tc6 = tt6 ./ (10 .^ D6[1:37,2])


Da1 = tt1 ./ tc1
Da2 = tt2 ./ tc2
Da3 = tt3 ./ tc3
Da4 = tt4 ./ tc4
Da5 = tt5 ./ tc5
Da6 = tt1 ./ tc6


plot(u1[1:37,1],log.(10,tc1[1:37]),color = :red3,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(u2[1:37,1],log.(10,tc2[1:37]),color = :blue3,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
plot!(u3[1:37,1],log.(10,tc3[1:37]),color = :green3,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:20,xlims=[0,20])
plot!(u4[1:37,1],log.(10,tc4[1:37]),color = :black,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
plot!(u5[1:37,1],log.(10,tc5[1:37]),color = :purple,linestyle =:dash,linewidth=2,markershape = :circle,markersize = :2,
labels=false,grid=false,thickness_scaling=1.5)
plot!(u6[1:37,1],log.(10,tc6[1:37]),color = :yellow4,linestyle =:dash,linewidth=2,markershape = :square,markersize = :1.5,
labels=false,grid=false,thickness_scaling=1.5)
# yticks=0.3:0.1:1,ylims=[0.3,1]
=#


#=
u1 = readdlm("Post_Process_Scalar/Data/Da_V1.txt",skipstart=1)
u2 = readdlm("Post_Process_Scalar/Data/Da_H1.txt",skipstart=1)
u3 = readdlm("Post_Process_Scalar/Data/Da_M1_prime.txt",skipstart=1)
u4 = readdlm("Post_Process_Scalar/Data/Da_H1_prime.txt",skipstart=1)
u5 = readdlm("Post_Process_Scalar/Data/Da_M1.txt",skipstart=1)
u6 = readdlm("Post_Process_Scalar/Data/Da_L1.txt",skipstart=1)

plot(u1[:,1],u1[:,2],color = :red3,linestyle =:dash,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(u2[:,1],u2[:,2],color = :blue3,linestyle =:dot,linewidth=2.5,
labels="H-L",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
plot!(u3[:,1],u3[:,2],color = :green3,linestyle =:dashdot,linewidth=2,
labels="H'-L",grid=false,thickness_scaling=1.5,xticks=0:5:21,xlims=[0,21])
plot!(u4[:,1],u4[:,2],color = :black,linestyle =:solid,linewidth=2,
labels="M'-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(u5[:,1],u5[:,2],color = :purple,linestyle =:dash,linewidth=2,markershape = :circle,markersize = :2,
labels="M-L",grid=false,thickness_scaling=1.5)
plot!(u6[:,1],u6[:,2],color = :yellow4,linestyle =:dash,linewidth=2,markershape = :square,markersize = :1.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#

#=
u1 = readdlm("Post_Process_Scalar/Data/Da_V2.txt",skipstart=1)
u2 = readdlm("Post_Process_Scalar/Data/Da_H2.txt",skipstart=1)
u3 = readdlm("Post_Process_Scalar/Data/Da_H2_prime.txt",skipstart=1)
u4 = readdlm("Post_Process_Scalar/Data/Da_M2_prime.txt",skipstart=1)
u5 = readdlm("Post_Process_Scalar/Data/Da_M2.txt",skipstart=1)
u6 = readdlm("Post_Process_Scalar/Data/Da_L2.txt",skipstart=1)

plot(u1[:,1],u1[:,2],color = :red3,linestyle =:dash,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(u2[:,1],u2[:,2],color = :blue3,linestyle =:dot,linewidth=2.5,
labels="H-L",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
plot!(u3[:,1],u3[:,2],color = :green3,linestyle =:dashdot,linewidth=2,
labels="H'-L",grid=false,thickness_scaling=1.5,xticks=0:5:21,xlims=[0,21])
plot!(u4[:,1],u4[:,2],color = :black,linestyle =:solid,linewidth=2,
labels="M'-L",grid=false,thickness_scaling=1.5)
plot!(u5[:,1],u5[:,2],color = :purple,linestyle =:dash,linewidth=2,markershape = :circle,markersize = :2,
labels="M-L",grid=false,thickness_scaling=1.5)
plot!(u6[:,1],u6[:,2],color = :yellow4,linestyle =:dash,linewidth=2,markershape = :square,markersize = :1.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#


# For Scalar forcing
supersat1 = readdlm("PR_DNS_Scalar/Case_L_2/supersat")
temperature1 = readdlm("PR_DNS_Scalar/Case_L_2/temperature")
RN1 = readdlm("PR_DNS_Scalar/Case_L_2/RN", skipstart=2)
transition1 = readdlm("PR_DNS_Scalar/Case_L_2/transition", skipstart=2)

# Parameter values
Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 0.0000216
sigma = 0.072
M_w = 0.018
R = 8.314
kappa = 0.61


# Calculate tao_microphysics

tao_c = zeros(80)
tao_turb = zeros(80)
Da = zeros(80)
rr = zeros(80)

r0 = RN1[1,3]
for i=1:80
    global T,r,S,G,A,Sk,tao_c,tao_turb,Da,rd,r0,dt,rr

    T = temperature1[i,2]
    r = RN1[i,3]
    S = supersat1[i,2]

    p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))

    G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
    G = 1/G

    A = (2*sigma*M_w)/(R*T*rho)
    rd = 1.2e-6
    Sk = (A/r) - kappa*(rd/r)^3
    tao_c[i] = r^2/abs(G*(S-Sk))
    tao_turb[i] = ((0.512^2)/transition1[i,2])^(1/3)
    Da[i] = tao_turb[i]/tao_c[i]

    # dt = 0.5
    # rr[i] = sqrt(r0^2 + 2*G*(S-Sk)*dt)
    # r0 = rr[i]
end

# plot(supersat1[1:80,1],log.(10,tao_turb),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

plot(supersat1[1:80,1],Da,color = :red4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))



supersat2 = readdlm("PR_DNS_Scalar/Case_L_L_2/supersat")
temperature2 = readdlm("PR_DNS_Scalar/Case_L_L_2/temperature")
RN2 = readdlm("PR_DNS_Scalar/Case_L_L_2/RN", skipstart=2)
transition2 = readdlm("PR_DNS_Scalar/Case_L_L_2/transition", skipstart=2)

# Calculate tao_microphysics

Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 0.0000216

tao_c = zeros(80)
tao_turb = zeros(80)
Da = zeros(80)

for i=1:80
    global T,r,S,G,A,Sk,tao_c,tao_turb,Da,rd

    T = temperature2[i,2]
    r = RN2[i,3]
    S = supersat2[i,2]

    p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))

    G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
    G = 1/G
    A = (2*sigma*M_w)/(R*T*rho)
    rd = 1.2e-6
    Sk = (A/r) - kappa*(rd/r)^3
    tao_c[i] = r^2/abs(G*(S-Sk))
    tao_turb[i] = ((0.512^2)/transition2[i,2])^(1/3)
    Da[i] = tao_turb[i]/tao_c[i]
end
Da[31] = (Da[30] + Da[32])/2

plot!(supersat2[1:80,1],Da,color = :blue4,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40,legend=false)


supersat3 = readdlm("PR_DNS_Scalar/Case_L_H_2/supersat")
temperature3 = readdlm("PR_DNS_Scalar/Case_L_H_2/temperature")
RN3 = readdlm("PR_DNS_Scalar/Case_L_H_2/RN", skipstart=2)
transition3 = readdlm("PR_DNS_Scalar/Case_L_H_2/transition", skipstart=2)

# Calculate tao_microphysics

Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 0.0000216

tao_c = zeros(80)
tao_turb = zeros(80)
Da = zeros(80)

for i=1:80
    global T,r,S,G,A,Sk,tao_c,tao_turb,Da,rd

    T = temperature3[i,2]
    r = RN3[i,3]
    S = supersat3[i,2]

    p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))

    G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
    G = 1/G
    A = (2*sigma*M_w)/(R*T*rho)
    rd = 1.2e-6
    Sk = (A/r) - kappa*(rd/r)^3
    tao_c[i] = r^2/abs(G*(S-Sk))
    tao_turb[i] = ((0.512^2)/transition3[i,2])^(1/3)
    Da[i] = tao_turb[i]/tao_c[i]
end

plot!(supersat3[1:80,1],Da,color = :summer,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
# yticks=-1:1:7,ylims=[-1.5,7]


supersat4 = readdlm("PR_DNS_Scalar/Case_H_2/supersat")
temperature4 = readdlm("PR_DNS_Scalar/Case_H_2/temperature")
RN4 = readdlm("PR_DNS_Scalar/Case_H_2/RN", skipstart=2)
transition4 = readdlm("PR_DNS_Scalar/Case_H_2/transition", skipstart=2)

# Calculate tao_microphysics

Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 0.0000216

tao_c = zeros(80)
tao_turb = zeros(80)
Da = zeros(80)

for i=1:80
    global T,r,S,G,A,Sk,tao_c,tao_turb,Da,rd

    T = temperature4[i,2]
    r = RN4[i,3]
    S = supersat4[i,2]

    p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))

    G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
    G = 1/G
    A = (2*sigma*M_w)/(R*T*rho)
    rd = 1.2e-6
    Sk = (A/r) - kappa*(rd/r)^3
    tao_c[i] = r^2/abs(G*(S-Sk))
    tao_turb[i] = ((0.512^2)/transition4[i,2])^(1/3)
    Da[i] = tao_turb[i]/tao_c[i]
end
Da[19] = (Da[18] + Da[20])/2

plot!(supersat4[1:80,1],Da,color = :yellow4,linestyle =:solid,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5)



supersat5 = readdlm("PR_DNS_Scalar/Case_H_L_2/supersat")
temperature5 = readdlm("PR_DNS_Scalar/Case_H_L_2/temperature")
RN5 = readdlm("PR_DNS_Scalar/Case_H_L_2/RN", skipstart=2)
transition5 = readdlm("PR_DNS_Scalar/Case_H_L_2/transition", skipstart=2)

# Calculate tao_microphysics

Lh = 2.5*1e6
Rv = 461.5
k = 0.0238
rho = 1000.0
mu_v = 0.0000216

tao_c = zeros(80)
tao_turb = zeros(80)
Da = zeros(80)

for i=1:80
    global T,r,S,G,A,Sk,tao_c,tao_turb,Da,rd

    T = temperature5[i,2]
    r = RN5[i,3]
    S = supersat5[i,2]

    p_sat = 611.2*exp(17.67*(T-273.15)/(T-29.65))

    G = ((Lh*rho)/(k*T))*((Lh/(Rv*T))-1) + (rho*Rv*T)/(mu_v*p_sat)
    G = 1/G
    A = (2*sigma*M_w)/(R*T*rho)
    rd = 1.2e-6
    Sk = (A/r) - kappa*(rd/r)^3
    tao_c[i] = r^2/abs(G*(S-Sk))
    tao_turb[i] = ((0.512^2)/transition5[i,2])^(1/3)
    Da[i] = tao_turb[i]/tao_c[i]
end
# Da[30] = (Da[29] + Da[31])/2

plot!(supersat5[1:80,1],Da,color = :black,linestyle =:dash,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5)
scatter!(supersat5[1:5:80,1],Da[1:5:80],color = :black,markershape = :hex,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)



#=
# For mixed-phase

supersat_w1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/supersat_water")
supersat_i1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/supersat_ice")
temperature1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/temperature")
transition1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/transition",skipstart=1)
vapor1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/vapor")
mixed_count1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/mixed-phase-count",skipstart=1)
mixed_rad1 = readdlm("PR_DNS_Mixed/F3_Case_L_S/mixed-phase-radii",skipstart=1)

supersat_w2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/supersat_water")
supersat_i2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/supersat_ice")
temperature2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/temperature")
transition2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/transition",skipstart=1)
vapor2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/vapor")
mixed_count2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/mixed-phase-count",skipstart=1)
mixed_rad2 = readdlm("PR_DNS_Mixed/F3_Case_L_C/mixed-phase-radii",skipstart=1)

supersat_w3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/supersat_water")
supersat_i3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/supersat_ice")
temperature3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/temperature")
transition3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/transition",skipstart=1)
vapor3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/vapor")
mixed_count3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/mixed-phase-count",skipstart=1)
mixed_rad3 = readdlm("PR_DNS_Mixed/F3_Case_L_P/mixed-phase-radii",skipstart=1)

supersat_w4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/supersat_water")
supersat_i4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/supersat_ice")
temperature4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/temperature")
transition4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/transition",skipstart=1)
vapor4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/vapor")
mixed_count4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/mixed-phase-count",skipstart=1)
mixed_rad4 = readdlm("PR_DNS_Mixed/F3_Case_H_S/mixed-phase-radii",skipstart=1)

supersat_w5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/supersat_water")
supersat_i5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/supersat_ice")
temperature5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/temperature")
transition5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/transition",skipstart=1)
vapor5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/vapor")
mixed_count5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/mixed-phase-count",skipstart=1)
mixed_rad5 = readdlm("PR_DNS_Mixed/F3_Case_L_S_2/mixed-phase-radii",skipstart=1)

supersat_w6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/supersat_water")
supersat_i6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/supersat_ice")
temperature6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/temperature")
transition6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/transition",skipstart=1)
vapor6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/vapor")
mixed_count6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/mixed-phase-count",skipstart=1)
mixed_rad6 = readdlm("PR_DNS_Mixed/F3_Case_H_S_2/mixed-phase-radii",skipstart=1)


# Parameter values
Lh_w = 2.5*1e6
Lh_ice = 2.83*1e6
Rv = 461.5
k = 0.0238
rho_w = 1000.0
rho_ice = 917.0
mu_v = 0.000022
mu = 1.64e-5

# Calculate tao_microphysics
tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature1[i,2]
    r_w = mixed_rad1[i,2]
    r_ice = mixed_rad1[i,6]
    S_w = supersat_w1[i,2]
    S_ice = supersat_i1[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition1[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end

plot(supersat_i1[:,1],log.(10,Da_ice),color = :red4,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# a = [supersat_i1[9,1]; supersat_i1[17,1]; supersat_i1[33,1]; supersat_i1[65,1]; supersat_i1[97,1]]
# b = [Da_ice[9]; Da_ice[17]; Da_ice[33]; Da_ice[65]; Da_ice[97]]
# c = [Da_w[9]; Da_w[17]; Da_w[33]; Da_w[65]; Da_w[97]]
# plot(a,log.(10,b),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature2[i,2]
    r_w = mixed_rad2[i,2]
    r_ice = mixed_rad2[i,6]
    S_w = supersat_w2[i,2]
    S_ice = supersat_i2[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition2[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end

plot!(supersat_i2[:,1],log.(10,Da_ice),color = :blue4,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,legend=false)


# a = [supersat_i2[9,1]; supersat_i2[17,1]; supersat_i2[33,1]; supersat_i2[65,1]; supersat_i2[97,1]]
# b = [Da_ice[9]; Da_ice[17]; Da_ice[33]; Da_ice[65]; Da_ice[97]]
# c = [Da_w[9]; Da_w[17]; Da_w[33]; Da_w[65]; Da_w[97]]
# # plot!(a,log.(10,b),color = :blue4,linestyle =:dash,linewidth=2,
# # labels=false,grid=false,thickness_scaling=1.5,legend=false)


tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature3[i,2]
    r_w = mixed_rad3[i,2]
    r_ice = mixed_rad3[i,6]
    S_w = supersat_w3[i,2]
    S_ice = supersat_i3[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition3[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end

plot!(supersat_i3[:,1],log.(10,Da_ice),color = :green4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,xticks=2:8:64,xlims=[2,64])

# a = [supersat_i3[9,1]; supersat_i3[17,1]; supersat_i3[33,1]; supersat_i3[65,1]; supersat_i3[97,1]]
# b = [Da_ice[9]; Da_ice[17]; Da_ice[33]; Da_ice[65]; Da_ice[97]]
# c = [Da_w[9]; Da_w[17]; Da_w[33]; Da_w[65]; Da_w[97]]
# # plot!(a,log.(10,b),color = :green4,linestyle =:dot,linewidth=2.5,
# # labels=false,grid=false,thickness_scaling=1.5,xticks=0:8:52,xlims=[0,52])


tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature4[i,2]
    r_w = mixed_rad4[i,2]
    r_ice = mixed_rad4[i,6]
    S_w = supersat_w4[i,2]
    S_ice = supersat_i4[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition4[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end


plot!(supersat_i4[:,1],log.(10,Da_ice),color = :yellow4,linestyle =:dashdot,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5)

# a = [supersat_i4[9,1]; supersat_i4[17,1]; supersat_i4[33,1]; supersat_i4[65,1]; supersat_i4[97,1]]
# b = [Da_ice[9]; Da_ice[17]; Da_ice[33]; Da_ice[65]; Da_ice[97]]
# c = [Da_w[9]; Da_w[17]; Da_w[33]; Da_w[65]; Da_w[97]]
# plot!(a,log.(10,b),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)


tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature5[i,2]
    r_w = mixed_rad5[i,2]
    r_ice = mixed_rad5[i,6]
    S_w = supersat_w5[i,2]
    S_ice = supersat_i5[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition4[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end

plot!(supersat_i5[:,1],log.(10,Da_ice),color = :black,linestyle =:dash,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5)
scatter!(supersat_i5[1:5:end,1],log.(10,Da_ice[1:5:end]),color = :black,
labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)


tao_w = zeros(129)
tao_ice = zeros(129)
Da_w = zeros(129)
Da_ice = zeros(129)

for i=3:129
    global T, r_w, r_ice, S_w, S_ice, p_sat_w, p_sat_ice, G_w, G_ice, tao_w, tao_ice, tao_turb, Da_w, Da_ice

    T = temperature6[i,2]
    r_w = mixed_rad6[i,2]
    r_ice = mixed_rad6[i,6]
    S_w = supersat_w6[i,2]
    S_ice = supersat_i6[i,2]

    p_sat_w = 611.2*exp(17.67*(T-273.15)/(T-29.65))
    p_sat_ice = 611.2*exp(22.46*(T-273.15)/(T-0.53))

    G_w = ((Lh_w*rho_w)/(k*T))*((Lh_w/(Rv*T))-1) + (rho_w*Rv*T)/(mu_v*p_sat_w)
    G_w = 1/G_w
    G_ice = ((Lh_ice*rho_ice)/(k*T))*((Lh_ice/(Rv*T))-1) + (rho_ice*Rv*T)/(mu_v*p_sat_ice)
    G_ice = 1/G_ice
    tao_w[i] = r_w^2/abs(G_w*S_w)
    tao_ice[i] = r_ice^2/abs(G_ice*S_ice)
    l = 0.2
    tao_turb = ((l^2)/transition4[i,2])^(1/3)
    Da_w[i] = tao_turb/tao_w[i]
    Da_ice[i] = tao_turb/tao_ice[i]
end

plot!(supersat_i6[:,1],log.(10,Da_ice),color = :purple,linestyle =:dot,linewidth=2.5,
labels="false",grid=false,thickness_scaling=1.5)
scatter!(supersat_i6[3:5:end,1],log.(10,Da_ice[3:5:end]),color = :purple,
labels=false,grid=false,thickness_scaling=1.5,markershape = :square,markersize = :1.5)
=#