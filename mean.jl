using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
#=
RN1 = readdlm("PR_DNS_Scalar/Case_L_2/RN", skipstart=2)
supersat1 = readdlm("PR_DNS_Scalar/Case_L_1/supersat")
temperature1 = readdlm("PR_DNS_Scalar/Case_L_1/temperature")
transition1 = readdlm("PR_DNS_Scalar/Case_L_1/transition", skipstart=1)
vapor1 = readdlm("PR_DNS_Scalar/Case_L_1/vapor")
Cd1 = readdlm("PR_DNS_Scalar/Case_L_1/Cd", skipstart=1)
Stat_1 = readdlm("PR_DNS_Scalar/Case_L_1/more_statistics", skipstart=1)

RN2 = readdlm("PR_DNS_Scalar/Case_L_L_2/RN", skipstart=2)
supersat2 = readdlm("PR_DNS_Scalar/Case_L_L_1/supersat")
temperature2 = readdlm("PR_DNS_Scalar/Case_L_L_1/temperature")
transition2 = readdlm("PR_DNS_Scalar/Case_L_L_1/transition", skipstart=1)
vapor2 = readdlm("PR_DNS_Scalar/Case_L_L_1/vapor")
Cd2 = readdlm("PR_DNS_Scalar/Case_L_L_1/Cd", skipstart=1)
Stat_2 = readdlm("PR_DNS_Scalar/Case_L_L_1/more_statistics", skipstart=1)

RN3 = readdlm("PR_DNS_Scalar/Case_L_H_2/RN", skipstart=2)
supersat3 = readdlm("PR_DNS_Scalar/Case_L_H_1/supersat")
temperature3 = readdlm("PR_DNS_Scalar/Case_L_H_1/temperature")
transition3 = readdlm("PR_DNS_Scalar/Case_L_H_1/transition", skipstart=1)
vapor3 = readdlm("PR_DNS_Scalar/Case_L_H_1/vapor")
Cd3 = readdlm("PR_DNS_Scalar/Case_L_H_1/Cd", skipstart=1)
Stat_3 = readdlm("PR_DNS_Scalar/Case_L_H_1/more_statistics", skipstart=1)

RN4 = readdlm("PR_DNS_Scalar/Case_H_2/RN", skipstart=2)
supersat4 = readdlm("PR_DNS_Scalar/Case_H_1/supersat")
temperature4 = readdlm("PR_DNS_Scalar/Case_H_1/temperature")
transition4 = readdlm("PR_DNS_Scalar/Case_H_1/transition", skipstart=1)
vapor4 = readdlm("PR_DNS_Scalar/Case_H_1/vapor")
Cd4 = readdlm("PR_DNS_Scalar/Case_H_1/Cd", skipstart=1)
Stat_4 = readdlm("PR_DNS_Scalar/Case_H_1/more_statistics", skipstart=1)

RN5 = readdlm("PR_DNS_Scalar/Case_H_L_2/RN", skipstart=2)
supersat5 = readdlm("PR_DNS_Scalar/Case_H_L_1/supersat")
temperature5 = readdlm("PR_DNS_Scalar/Case_H_L_1/temperature")
transition5 = readdlm("PR_DNS_Scalar/Case_H_L_1/transition", skipstart=1)
vapor5 = readdlm("PR_DNS_Scalar/Case_H_L_1/vapor")
Cd5 = readdlm("PR_DNS_Scalar/Case_H_L_1/Cd", skipstart=1)
Stat_5 = readdlm("PR_DNS_Scalar/Case_H_L_1/more_statistics", skipstart=1)

RN6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/RN", skipstart=2)
supersat6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/supersat")
temperature6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/temperature")
transition6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/transition", skipstart=1)
vapor6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/vapor")
Cd6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/Cd", skipstart=2)
Stat_6 = readdlm("PR_DNS_Scalar/Case_L_ALR_1/more_statistics", skipstart=1)
=#

# plot(temperature1[:,1],temperature1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(temperature2[:,1],temperature2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="false",legend=false,xticks=0:5:40)
# plot!(temperature3[:,1],temperature3[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(temperature4[:,1],temperature4[:,2],color = :orange3,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(temperature5[:,1],temperature5[:,2],color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false)
# scatter!(temperature5[1:6:end,1],temperature5[1:6:end,2],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(temperature6[:,1],temperature6[:,2],color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(temperature6[1:5:end,1],temperature6[1:5:end,2],color = :purple,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=270.74:0.02:270.88,ylims=[270.74,270.89]


# plot(temperature1[:,1],log.(10,temperature1[:,3]),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(temperature2[:,1],log.(10,temperature2[:,3]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(temperature3[:,1],log.(10,temperature3[:,3]),color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(temperature4[:,1],log.(10,temperature4[:,3]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(temperature5[:,1],log.(10,temperature5[:,3]),color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(temperature5[1:3:end,1],log.(10,temperature5[1:3:end,3]),color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(temperature6[:,1],log.(10,temperature6[:,3]),color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(temperature6[1:5:end,1],log.(10,temperature6[1:5:end,3]),color = :purple,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=-10:2:-2,ylims=[-10,-2]


# plot(vapor1[:,1],vapor1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(vapor2[:,1],vapor2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="false",legend=false,xticks=0:5:40)
# plot!(vapor3[:,1],vapor3[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(vapor4[:,1],vapor4[:,2],color = :orange3,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(vapor5[:,1],vapor5[:,2],color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,yticks=3.87:0.02:3.97,ylims=[3.868,3.972])
# scatter!(vapor5[1:6:end,1],vapor5[1:6:end,2],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)

# yticks=3.87:0.02:3.97,ylims=[3.868,3.972]


# plot(vapor1[:,1],log.(10,vapor1[:,3]),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(vapor2[:,1],log.(10,vapor2[:,3]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(vapor3[:,1],log.(10,vapor3[:,3]),color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(vapor4[:,1],log.(10,vapor4[:,3]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(vapor5[:,1],log.(10,vapor5[:,3]),color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(vapor5[1:3:end,1],log.(10,vapor5[1:3:end,3]),color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# # yticks=-8:1.5:1,ylims=[-8,1]


# plot(supersat1[:,1],supersat1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],supersat2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="false",legend=false,xticks=0:5:40)
# plot!(supersat3[:,1],supersat3[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(supersat4[:,1],supersat4[:,2],color = :orange3,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(supersat5[:,1],supersat5[:,2],color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(supersat5[1:6:end,1],supersat5[1:6:end,2],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(supersat6[:,1],supersat6[:,2],color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(supersat6[1:6:end,1],supersat6[1:6:end,2],color = :purple,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=-0.5:0.1:0.0,ylims=[-0.5,0.01]


# plot(supersat1[:,1],log.(10,supersat1[:,3]),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],log.(10,supersat2[:,3]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(supersat3[:,1],log.(10,supersat3[:,3]),color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(supersat4[:,1],log.(10,supersat4[:,3]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(supersat5[:,1],log.(10,supersat5[:,3]),color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(supersat5[1:5:end,1],log.(10,supersat5[1:5:end,3]),color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(supersat6[:,1],log.(10,supersat6[:,3]),color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(supersat6[1:5:end,1],log.(10,supersat6[1:5:end,3]),color = :purple,markershape = :circle,markersize = :2,
# labels=false,grid=false,thickness_scaling=1.5)


# plot(Cd1[:,1],Cd1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(Cd2[:,1],Cd2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="false",legend=false,xticks=0:5:40)
# plot!(Cd3[:,1],Cd3[:,2],color = :summer,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(Cd4[:,1],Cd4[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)


# plot(RN1[:,1],RN1[:,3] .*1e6,color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[:,1],RN2[:,3] .*1e6,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(RN3[:,1],RN3[:,3] .*1e6,color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(RN4[:,1],RN4[:,3] .*1e6,color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(RN5[:,1],RN5[:,3] .*1e6,color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN5[1:5:end,1],RN5[1:5:end,3] .*1e6,color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(RN6[:,1],RN6[:,3] .*1e6,color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN6[1:5:end,1],RN6[1:5:end,3] .*1e6,color = :purple,markershape = :circle,markersize = :2,
# labels=false,grid=false,yticks=0:1:5,ylims=[0,5])
# yticks=3:2:15,ylims=[3,15.1]
# yticks=0:1:5,ylims=[0,5]

# plot(RN1[:,1],log.(10,RN1[:,4]),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[:,1],log.(10,RN2[:,4]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(RN3[:,1],log.(10,RN3[:,4]),color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(RN4[:,1],log.(10,RN4[:,4]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(RN5[:,1],log.(10,RN5[:,4]),color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN5[1:6:end,1],log.(10,RN5[1:6:end,4]),color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(RN6[:,1],log.(10,RN6[:,4]),color = :purple,linestyle =:dot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN6[1:6:end,1],log.(10,RN6[1:6:end,4]),color = :purple,markershape = :utriangle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=-7.5:0.5:-5.5,ylims=[-7.7,-5.5]

## Activation/Deactivation Fraction
# N = 16e6
# plot(Stat_1[:,1],(Stat_1[:,3]) .*1e-6,color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(Stat_2[:,1],(Stat_2[:,3]) .*1e-6,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(Stat_3[:,1],(Stat_3[:,3]) .*1e-6,color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(Stat_4[:,1],(Stat_4[:,3]) .*1e-6,color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(Stat_5[:,1],(Stat_5[:,3]) .*1e-6,color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# scatter!(Stat_5[1:6:end,1],(Stat_5[1:6:end,3]) .*1e-6,color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(Stat_6[:,1],(Stat_6[:,3]) .*1e-6,color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(Stat_6[1:6:end,1],(Stat_6[1:6:end,3]) .*1e-6,color = :purple,markershape = :circle,markersize = :2,
# labels=false)

# plot(Stat_1[40:end,1],(Stat_1[40:end,3]) .*1e-6,color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(Stat_2[40:end,1],(Stat_2[40:end,3]) .*1e-6,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(Stat_3[40:end,1],(Stat_3[40:end,3]) .*1e-6,color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(Stat_4[40:end,1],(Stat_4[40:end,3]) .*1e-6,color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(Stat_5[40:end,1],(Stat_5[40:end,3]) .*1e-6,color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# scatter!(Stat_5[40:6:end,1],(Stat_5[40:6:end,3]) .*1e-6,color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(Stat_6[40:end,1],(Stat_6[40:end,3]) .*1e-6,color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(Stat_6[40:6:end,1],(Stat_6[40:6:end,3]) .*1e-6,color = :purple,markershape = :circle,markersize = :2,
# labels=false)


# # <r'S'>
# plot(Stat_1[:,1],Stat_1[:,6],color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(Stat_2[:,1],Stat_2[:,6],color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(Stat_3[:,1],Stat_3[:,6],color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false)
# plot!(Stat_4[:,1],Stat_4[:,6],color = :orange3,linestyle =:solid,linewidth=2,
# labels="L-LF",grid=false,thickness_scaling=1.5,legend=false)
# plot!(Stat_5[:,1],Stat_5[:,6],color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# scatter!(Stat_5[1:6:end,1],Stat_5[1:6:end,6],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(Stat_6[:,1],Stat_6[:,6],color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false)
# scatter!(Stat_6[1:6:end,1],Stat_6[1:6:end,6],color = :purple,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=-0.02:0.02:0.06,ylims=[-0.02,0.06]


# T1 = zeros(81)
# T1[1] = temperature1[1,2]
# T2 = zeros(81)
# T2[1] = temperature2[1,2]

# m = (2.5e6)/1005
# for i=1:80
#     T1[i+1] = m*(Cd1[i+1,2]+Cd1[i,2])/2 *(Cd1[i+1,1]-Cd1[i,1]) + T1[i]
#     T2[i+1] = m*(Cd2[i+1,2]+Cd2[i,2])/2 *(Cd2[i+1,1]-Cd2[i,1]) + T2[i]
# end

# Conv_Diff = temperature1[:,2] .- T1
# temperature2[:,2] = temperature2[:,2] .- Conv_Diff

# plot(Cd1[:,1],T1[:],color = :bwr,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40,yticks=269.75:0.25:271,ylims=[269.75,271.1])
# plot!(Cd1[:,1],temperature1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(Cd2[:,1],T2[:],color = :bwr,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40,yticks=269.25:0.25:271,ylims=[269.25,271.1])
# plot!(Cd2[:,1],temperature2[:,2],color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# Relative dipsersions
# c1 = sqrt.(temperature1[:,3])./temperature1[:,2]
# c2 = sqrt.(temperature2[:,3])./temperature2[:,2]
# c3 = sqrt.(temperature3[:,3])./temperature3[:,2]
# c4 = sqrt.(temperature4[:,3])./temperature4[:,2]
# c5 = sqrt.(temperature5[:,3])./temperature5[:,2]

# c1 = sqrt.(vapor1[:,3])./vapor1[:,2]
# c2 = sqrt.(vapor2[:,3])./vapor2[:,2]
# c3 = sqrt.(vapor3[:,3])./vapor3[:,2]
# c4 = sqrt.(vapor4[:,3])./vapor4[:,2]
# c5 = sqrt.(vapor5[:,3])./vapor5[:,2]

# c1 = sqrt.(supersat1[:,3])./abs.(supersat1[:,2])
# c2 = sqrt.(supersat2[:,3])./abs.(supersat2[:,2])
# c3 = sqrt.(supersat3[:,3])./abs.(supersat3[:,2])
# c4 = sqrt.(supersat4[:,3])./abs.(supersat4[:,2])
# c5 = sqrt.(supersat5[:,3])./abs.(supersat5[:,2])

# plot(supersat1[:,1],log.(10,c1),color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],log.(10,c2),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(supersat3[:,1],log.(10,c3),color = :summer,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(supersat4[:,1],log.(10,c4),color = :yellow4,linestyle =:solid,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5,yticks=-3:1:2,ylims=[-3,2])
# plot!(supersat5[:,1],log.(10,c5),color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,yticks=-7:1:-3,ylims=[-7,-2.6])
# scatter!(supersat5[1:3:end,1],log.(10,c5[1:3:end,1]),color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# # yticks=-2.5:1:1.5,ylims=[-2.5,1.5]

# c1 = RN1[:,4]./RN1[:,3]
# c2 = RN2[:,4]./RN2[:,3]
# c3 = RN3[:,4]./RN3[:,3]
# c4 = RN4[:,4]./RN4[:,3]
# c5 = RN5[:,4]./RN5[:,3]
# c6 = RN6[:,4]./RN6[:,3]

# c1 = RN1[:,4] .*1e6
# c2 = RN2[:,4] .*1e6
# c3 = RN3[:,4] .*1e6
# c4 = RN4[:,4] .*1e6
# c5 = RN5[:,4] .*1e6

# plot(RN1[:,1],c1,color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[:,1],c2,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(RN3[:,1],c3,color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40)
# plot!(RN4[:,1],c4,color = :orange3,linestyle =:solid,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(RN5[:,1],c5,color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN5[1:6:end,1],c5[1:6:end],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# plot!(RN6[:,1],c6,color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(RN6[1:6:end,1],c6[1:6:end],color = :purple,markershape = :circle,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5)
# yticks=0:1:5,ylims=[0,5]

# c1 = sqrt.(temperature1[:,3])./temperature1[:,2]
# c2 = sqrt.(temperature2[:,3])./temperature2[:,2]
# c3 = sqrt.(temperature3[:,3])./temperature3[:,2]
# c4 = sqrt.(temperature4[:,3])./temperature4[:,2]
# c5 = sqrt.(temperature5[:,3])./temperature5[:,2]
# c6 = sqrt.(temperature6[:,3])./temperature6[:,2]

# plot(transition1[:,1],c1,color = :red4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(transition2[:,1],c2,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:40,xlims=[0,40])
# plot!(transition3[:,1],c3,color = :green4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(transition4[:,1],c4,color = :orange3,linestyle =:solid,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(transition5[:,1],c5,color = :black,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(transition5[1:6:end,1],c5[1:6:end],color = :black,markershape = :square,markersize = :2,
# labels="Critical",grid=false,thickness_scaling=1.5,legend=false)
# plot!(transition6[:,1],c6,color = :purple,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(transition6[1:6:end,1],c6[1:6:end],color = :purple,markershape = :circle,markersize = :2,
# labels=false,grid=false,thickness_scaling=1.5)


#=
plot(transition1[:,1],transition1[:,2],color = :brown,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(transition2[:,1],transition2[:,2],color = :summer,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
plot!(transition3[:,1],transition3[:,2],color = :purple,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
plot!(transition4[:,1],transition4[:,2],color = :nuuk,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
=#

#=
val1 = readdlm("Post_Process_Scalar/Data/L.txt",skipstart=1)
val2 = readdlm("Post_Process_Scalar/Data/M.txt",skipstart=1)
val3 = readdlm("Post_Process_Scalar/Data/H.txt",skipstart=1)
val4 = readdlm("Post_Process_Scalar/Data/V.txt",skipstart=1)

plot(val1[:,1],val1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="0.004,4.12",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(val2[:,1],val2[:,2],color = :black,linestyle =:dashdot,linewidth=2,
labels="0.004,8.12",grid=false,thickness_scaling=1.5)
plot!(val3[:,1],val3[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="0.016,4.12",grid=false,thickness_scaling=1.5)
plot!(val4[:,1],val4[:,2],color = :green4,linestyle =:solid,linewidth=2,
labels="0.016,8.12",grid=false,thickness_scaling=1.5)
=#

#=
A1 = readdlm("Post_Process_Scalar/Data/v_V_2.txt",skipstart=1)
A1[15:end,2] .= 0
A2 = readdlm("Post_Process_Scalar/Data/v_H_2.txt",skipstart=1)
A2[5:end,2] .= 0
new_row = [24.0 0]
new_row2 = [26.0 0]
vcat(A2, new_row)
vcat(A2, new_row2)
A3 = readdlm("Post_Process_Scalar/Data/v_M_2.txt",skipstart=1)
A3[9:end,2] .= 0
vcat(A3, new_row)
vcat(A3, new_row2)
A4 = readdlm("Post_Process_Scalar/Data/v_L_2.txt",skipstart=1)
A4[9:end,2] .= 0
vcat(A4, new_row)
vcat(A4, new_row2)


plot(A1[1:end-2,1] .-A1[1,1],A1[1:end-2,2],color = :red3,linestyle =:dash,linewidth=2,markershape=:circle,markersize=:2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(A2[1:end,1] .-A2[1,1],A2[1:end,2],color = :blue3,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
plot!(A3[1:end,1] .-A3[1,1],A3[1:end,2],color = :green3,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,xticks=0:5:20,xlims=[-0.5,20])
plot!(A4[1:end,1] .-A4[1,1],A4[1:end,2],color = :yellow4,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
# # ylims=[0.002,0.01],yticks=0.002:0.001:0.01,
=#

#=
transition5 = readdlm("PR_DNS_Cloud/F3_CaseM1_128/transition", skipstart=2)
transition6 = readdlm("PR_DNS_Cloud/F3_CaseH1_128/transition", skipstart=2)

plot(transition5[1:30,1],transition5[1:30,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="0.004,3.91",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(transition6[1:30,1],transition6[1:30,2],color = :black,linestyle =:dashdot,linewidth=2,
labels="0.004,8.11",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))

A1 = readdlm("Post_Process_Scalar/Data/M_eps.txt",skipstart=1)
A2 = readdlm("Post_Process_Scalar/Data/H_eps.txt",skipstart=1)

plot!(A1[:,1],A1[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="0.016,3.91",grid=false,thickness_scaling=1.5,ylims=[0.0,0.005],xticks=0.5:3:16,xlims=[0,15])
plot!(A2[:,1],A2[:,2],color = :green4,linestyle =:solid,linewidth=2,
labels="0.016,8.11",grid=false,thickness_scaling=1.5)
=#

#=
A1 = readdlm("Post_Process_Scalar/Data/a.txt",skipstart=1)
A2 = readdlm("Post_Process_Scalar/Data/c.txt",skipstart=1)
A3 = readdlm("Post_Process_Scalar/Data/b.txt",skipstart=1)
A4 = readdlm("Post_Process_Scalar/Data/d.txt",skipstart=1)

plot(A1[1:end-7,1],A1[1:end-7,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(A2[1:end-6,1],A2[1:end-6,2],color = :black,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
plot!(A3[1:end-7,1],A3[1:end-7,2],color = :blue4,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,ylims=[-12,-2],yticks=-12:2:-2,xlims=[-2,0.5])
plot!(A4[1:end-6,1],A4[1:end-6,2],color = :green4,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
=#

#=
A1 = readdlm("Post_Process_Scalar/Data/urms_V1.txt",skipstart=1)
A2 = readdlm("Post_Process_Scalar/Data/urms_H1.txt",skipstart=1)
A3 = readdlm("Post_Process_Scalar/Data/urms_M2.txt",skipstart=1)
A4 = readdlm("Post_Process_Scalar/Data/urms_L2.txt",skipstart=1)

A2[1,2] = A1[1,2]
A3[1,2] = A1[1,2]
A4[1,2] = A1[1,2]

plot(A1[:,1],A1[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(A2[:,1],A2[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,ylims=[0.04,0.22],yticks=0.04:0.02:0.22,xlims=[0,24])
plot!(A3[:,1],A3[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="M-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(A4[:,1],A4[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#


#=
A1 = readdlm("PR_DNS_Cloud/F3_CaseV1_128/PDF-zvel/zvel-10.00",skipstart=1)
A2 = readdlm("PR_DNS_Cloud/F3_CaseH1_128/PDF-zvel/zvel-10.00",skipstart=1)
A3 = readdlm("PR_DNS_Cloud/F3_CaseM1_128/PDF-zvel/zvel-10.00",skipstart=1)
A4 = readdlm("PR_DNS_Cloud/F3_CaseL1_128/PDF-zvel/zvel-10.00",skipstart=1)

plot(A1[:,1],A1[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(A2[:,1],A2[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5)
plot!(A3[:,1],A3[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="M-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(A4[:,1],A4[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#

#=
A1 = readdlm("PR_DNS_Cloud/F3_CaseV1_128/PDF-temperature/temperature-10.00",skipstart=1)
A2 = readdlm("PR_DNS_Cloud/F3_CaseH1_128/PDF-temperature/temperature-10.00",skipstart=1)
A3 = readdlm("PR_DNS_Cloud/F3_CaseM1_128/PDF-temperature/temperature-10.00",skipstart=1)
A4 = readdlm("PR_DNS_Cloud/F3_CaseL1_128/PDF-temperature/temperature-10.00",skipstart=1)

M1 = readdlm("PR_DNS_Cloud/F3_CaseV1_128/temperature",skipstart=1)
mean1 = M1[20,2]
M2 = readdlm("PR_DNS_Cloud/F3_CaseH1_128/temperature",skipstart=1)
mean2 = M2[20,2]
M3 = readdlm("PR_DNS_Cloud/F3_CaseM1_128/temperature",skipstart=1)
mean3 = M3[20,2]
M4 = readdlm("PR_DNS_Cloud/F3_CaseL1_128/temperature",skipstart=1)
mean4 = M4[20,2]

plot(A1[:,1].- mean1,A1[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(A2[:,1].- mean2,A2[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,ylims=[0,800],yticks=0:200:800,xlims=[-0.015,0.015],xticks=-0.015:0.0075:0.015)
plot!(A3[:,1].- mean3,A3[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="M-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(A4[:,1].- mean4,A4[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#
#=
A1 = readdlm("PR_DNS_Cloud/F3_CaseV2_128/PDF-vapor/vapor-10.00",skipstart=1)
A2 = readdlm("PR_DNS_Cloud/F3_CaseH2_128/PDF-vapor/vapor-10.00",skipstart=1)
A3 = readdlm("PR_DNS_Cloud/F3_CaseM2_128/PDF-vapor/vapor-10.00",skipstart=1)
A4 = readdlm("PR_DNS_Cloud/F3_CaseL2_128/PDF-vapor/vapor-10.00",skipstart=1)

M1 = readdlm("PR_DNS_Cloud/F3_CaseV2_128/vapor",skipstart=1)
mean1 = M1[20,2]
M2 = readdlm("PR_DNS_Cloud/F3_CaseH2_128/vapor",skipstart=1)
mean2 = M2[20,2]
M3 = readdlm("PR_DNS_Cloud/F3_CaseM2_128/vapor",skipstart=1)
mean3 = M3[20,2]
M4 = readdlm("PR_DNS_Cloud/F3_CaseL2_128/vapor",skipstart=1)
mean4 = M4[20,2]

# ylims=[0,700],yticks=0:100:700,xlims=[-0.015,0.015],xticks=-0.015:0.0075:0.015
plot(A1[:,1].- mean1,A1[:,2],color = :green4,linestyle =:dashdot,linewidth=2,
labels="V-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
plot!(A2[:,1].- mean2,A2[:,2],color = :yellow4,linestyle =:solid,linewidth=2,ylims=[0,10],yticks=0:2:10,xlims=[-1,0.6],xticks=-1:0.4:0.6,
labels="H-L",grid=false,thickness_scaling=1.5)
plot!(A3[:,1].- mean3,A3[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="M-L",grid=false,thickness_scaling=1.5,legend=false)
plot!(A4[:,1].- mean4,A4[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
labels="L-L",grid=false,thickness_scaling=1.5)
=#
#=
RN11 = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F/RN", skipstart=2)
RN22 = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F0/RN", skipstart=2)

supersat11 = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F/supersat")
supersat22 = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F0/supersat")

# plot(RN11[1:end-60,1],log.(10,RN11[1:end-60,4].^2),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN22[:,1],log.(10,RN22[:,4].^2),color = :summer,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)

plot(supersat11[1:end-60,1],log.(10,supersat11[1:end-60,3]),color = :blue4,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(supersat22[:,1],log.(10,supersat22[:,3]),color = :summer,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,legend=false)
=#

# RN11 = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F/RN", skipstart=2)
# plot(r,R,color = :orange,linestyle =:dash,linewidth=2,markershape = :circle,markersize = :2,
# labels="R_u'v', L",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"),tickfont=font(9,"Helvetica Bold"))

#=
## FOR MIXED-PHASE CLOUD STATISTICS

supersat_w1 = readdlm("PR_DNS_Mixed/Case_L_S_L/supersat_water")
supersat_i1 = readdlm("PR_DNS_Mixed/Case_L_S_L/supersat_ice")
temperature1 = readdlm("PR_DNS_Mixed/Case_L_S_L/temperature")
transition1 = readdlm("PR_DNS_Mixed/Case_L_S_L/transition",skipstart=1)
vapor1 = readdlm("PR_DNS_Mixed/Case_L_S_L/vapor")
mixed_count1 = readdlm("PR_DNS_Mixed/Case_L_S_L/mixed-phase-count",skipstart=1)
mixed_rad1 = readdlm("PR_DNS_Mixed/Case_L_S_L/mixed-phase-radii",skipstart=1)

supersat_w2 = readdlm("PR_DNS_Mixed/Case_L_C_L/supersat_water")
supersat_i2 = readdlm("PR_DNS_Mixed/Case_L_C_L/supersat_ice")
temperature2 = readdlm("PR_DNS_Mixed/Case_L_C_L/temperature")
transition2 = readdlm("PR_DNS_Mixed/Case_L_C_L/transition",skipstart=1)
vapor2 = readdlm("PR_DNS_Mixed/Case_L_C_L/vapor")
mixed_count2 = readdlm("PR_DNS_Mixed/Case_L_C_L/mixed-phase-count",skipstart=1)
mixed_rad2 = readdlm("PR_DNS_Mixed/Case_L_C_L/mixed-phase-radii",skipstart=1)

supersat_w3 = readdlm("PR_DNS_Mixed/Case_L_P_L/supersat_water")
supersat_i3 = readdlm("PR_DNS_Mixed/Case_L_P_L/supersat_ice")
temperature3 = readdlm("PR_DNS_Mixed/Case_L_P_L/temperature")
transition3 = readdlm("PR_DNS_Mixed/Case_L_P_L/transition",skipstart=1)
vapor3 = readdlm("PR_DNS_Mixed/Case_L_P_L/vapor")
mixed_count3 = readdlm("PR_DNS_Mixed/Case_L_P_L/mixed-phase-count",skipstart=1)
mixed_rad3 = readdlm("PR_DNS_Mixed/Case_L_P_L/mixed-phase-radii",skipstart=1)

supersat_w4 = readdlm("PR_DNS_Mixed/Case_H_S_L/supersat_water")
supersat_i4 = readdlm("PR_DNS_Mixed/Case_H_S_L/supersat_ice")
temperature4 = readdlm("PR_DNS_Mixed/Case_H_S_L/temperature")
transition4 = readdlm("PR_DNS_Mixed/Case_H_S_L/transition",skipstart=1)
vapor4 = readdlm("PR_DNS_Mixed/Case_H_S_L/vapor")
mixed_count4 = readdlm("PR_DNS_Mixed/Case_H_S_L/mixed-phase-count",skipstart=1)
mixed_rad4 = readdlm("PR_DNS_Mixed/Case_H_S_L/mixed-phase-radii",skipstart=1)

supersat_w5 = readdlm("PR_DNS_Mixed/Case_H_S_H/supersat_water")
supersat_i5 = readdlm("PR_DNS_Mixed/Case_H_S_H/supersat_ice")
temperature5 = readdlm("PR_DNS_Mixed/Case_H_S_H/temperature")
transition5 = readdlm("PR_DNS_Mixed/Case_H_S_H/transition",skipstart=1)
vapor5 = readdlm("PR_DNS_Mixed/Case_H_S_H/vapor")
mixed_count5 = readdlm("PR_DNS_Mixed/Case_H_S_H/mixed-phase-count",skipstart=1)
mixed_rad5 = readdlm("PR_DNS_Mixed/Case_H_S_H/mixed-phase-radii",skipstart=1)

# supersat_w6 = readdlm("PR_DNS_Mixed/Case_HHS/supersat_water")
# supersat_i6 = readdlm("PR_DNS_Mixed/Case_HHS/supersat_ice")
# temperature6 = readdlm("PR_DNS_Mixed/Case_HHS/temperature")
# transition6 = readdlm("PR_DNS_Mixed/Case_HHS/transition",skipstart=1)
# vapor6 = readdlm("PR_DNS_Mixed/Case_HHS/vapor")
# mixed_count6 = readdlm("PR_DNS_Mixed/Case_HHS/mixed-phase-count",skipstart=1)
# mixed_rad6 = readdlm("PR_DNS_Mixed/Case_HHS/mixed-phase-radii",skipstart=1)
=#

# yticks=259.35:0.05:259.57,ylims=[259.35,259.57]
# plot(temperature1[:,1],temperature1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels="L-H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(temperature2[:,1],temperature2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="H-H",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(temperature3[:,1],temperature3[:,2],color = :green4,linestyle =:dot,linewidth=2.5,
# labels="L-L",grid=false,thickness_scaling=1.5,xticks=0:8:52,xlims=[0,52])
# plot!(temperature4[:,1],temperature4[:,2],color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="H-L",grid=false,thickness_scaling=1.5,legend=false)
# plot!(temperature5[:,1],temperature5[:,2],color = :black,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(temperature5[1:3:end,1],temperature5[1:3:end,2],color = :black,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(temperature6[:,1],temperature6[:,2],color = :purple,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(temperature6[1:3:end,1],temperature6[1:3:end,2],color = :purple,markershape = :square,markersize = :1.5,
# labels="false",grid=false,thickness_scaling=1.5)


# # yticks=2.200:0.025:2.325,ylims=[2.200,2.325]
# plot(vapor1[:,1],vapor1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(vapor2[:,1],vapor2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5,legend=false)
# plot!(vapor3[:,1],vapor3[:,2],color = :green4,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(vapor4[:,1],vapor4[:,2],color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(vapor5[:,1],vapor5[:,2],color = :black,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(vapor5[1:5:end,1],vapor5[1:5:end,2],color = :black,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(vapor6[:,1],vapor6[:,2],color = :purple,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(vapor6[1:3:end,1],vapor6[1:3:end,2],color = :purple,markershape = :square,markersize = :1.5,
# labels="false",grid=false,thickness_scaling=1.5)


# plot(supersat_w1[:,1],supersat_w1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat_w2[:,1],supersat_w2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(supersat_w3[:,1],supersat_w3[:,2],color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(supersat_w4[:,1],supersat_w4[:,2],color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,yticks=-0.09:0.02:-0.005,ylims=[-0.09,-0.005])
# plot!(supersat_w5[:,1],supersat_w5[:,2],color = :orange3,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(supersat_w5[1:5:end,1],supersat_w5[1:5:end,2],color = :orange3,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)

### AIAA 2027 #Start
# plot(supersat_w1[:,1],supersat_w1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat_w4[:,1],supersat_w4[:,2],color = :green4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,yticks=-0.09:0.02:-0.005,ylims=[-0.09,-0.005])
# plot!(supersat_w5[:,1],supersat_w5[:,2],color = :black,linestyle =:dot,linewidth=2.5,
# legend=false,grid=false,thickness_scaling=1.5,xticks=0:10:60,xlims=[0,60])

# plot(supersat_i1[:,1],supersat_i1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat_i4[:,1],supersat_i4[:,2],color = :green4,linestyle =:dash,linewidth=2,
# labels=false,grid=false)
# plot!(supersat_i5[:,1],supersat_i5[:,2],color = :black,linestyle =:dot,linewidth=2.5,
# legend=false,grid=false,thickness_scaling=1.5,xticks=0:10:60,xlims=[0,60])

# plot(mixed_rad1[:,1],mixed_rad1[:,2] .*1e6,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_rad4[:,1],mixed_rad4[:,2] .*1e6,color = :green4,linestyle =:dash,linewidth=2,
# labels=false,grid=false)
# plot!(mixed_rad5[:,1],mixed_rad5[:,2] .*1e6,color = :black,linestyle =:dot,linewidth=2.5,
# legend=false,grid=false,thickness_scaling=1.5,xticks=0:10:60,xlims=[0,60])

# plot(mixed_rad1[:,1],mixed_rad1[:,4] .*1e6,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_rad4[:,1],mixed_rad4[:,4] .*1e6,color = :green4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,ylims=[0,20])
# plot!(mixed_rad5[:,1],mixed_rad5[:,4] .*1e6,color = :black,linestyle =:dot,linewidth=2.5,
# legend=false,grid=false,thickness_scaling=1.5,xticks=0:10:60,xlims=[0,60])

### AIAA 2027 #End


# v1 = zeros(241); v2 = zeros(241)
# plot!(twinx(),v1,v2,color = :red4,linestyle =:solid,linewidth=2,yticks=0:0.02:0.085,ylims=[-0.007,0.085],
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# plot(supersat_i1[:,1],supersat_i1[:,2],color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat_i2[:,1],supersat_i2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(supersat_i3[:,1],supersat_i3[:,2],color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(supersat_i4[:,1],supersat_i4[:,2],color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,)
# plot!(supersat_i5[:,1],supersat_i5[:,2],color = :orange3,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(supersat_i5[1:5:end,1],supersat_i5[1:5:end,2],color = :orange3,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)

# plot(supersat_i1[:,1],log.(10,supersat_i1[:,3]),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat_i2[:,1],log.(10,supersat_i2[:,3]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(supersat_i3[:,1],log.(10,supersat_i3[:,3]),color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(supersat_i4[:,1],log.(10,supersat_i4[:,3]),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,)
# plot!(supersat_i5[:,1],log.(10,supersat_i5[:,3]),color = :orange3,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(supersat_i5[1:5:end,1],log.(10,supersat_i5[1:5:end,3]),color = :orange3,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)



## Droplet radii
# plot(mixed_rad1[:,1],mixed_rad1[:,2] .*1e6,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_rad2[:,1],mixed_rad2[:,2] .*1e6,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_rad3[:,1],mixed_rad3[:,2] .*1e6,color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_rad4[:,1],mixed_rad4[:,2] .*1e6,color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(mixed_rad5[:,1],mixed_rad5[:,2] .*1e6,color = :orange3,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(mixed_rad5[1:5:end,1],mixed_rad5[1:5:end,2] .*1e6,color = :orange3,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(mixed_rad6[:,1],mixed_rad6[:,2] .*1e6,color = :purple,linestyle =:dot,linewidth=2.5,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(mixed_rad6[1:10:end,1],mixed_rad6[1:10:end,2] .*1e6,color = :purple,markershape = :star4,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5)
# yticks=1:1:8,ylims=[1,8]

## Ice radii
# plot(mixed_rad1[:,1],mixed_rad1[:,4] .*1e6,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_rad2[:,1],mixed_rad2[:,4] .*1e6,color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_rad3[:,1],mixed_rad3[:,4] .*1e6,color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_rad4[:,1],mixed_rad4[:,4] .*1e6,color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(mixed_rad5[:,1],mixed_rad5[:,4] .*1e6,color = :orange3,linestyle =:dash,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# scatter!(mixed_rad5[1:5:end,1],mixed_rad5[1:5:end,4] .*1e6,color = :orange3,markershape = :circle,markersize = :2,
# labels="false",grid=false,thickness_scaling=1.5,yticks=0:5:35,ylims=[0,35])


## Ice mass fraction
# plot(mixed_count1[:,1],mixed_count1[:,5] ./(mixed_count1[:,4] .+ mixed_count1[:,5]),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[:,1],mixed_count2[:,5] ./(mixed_count2[:,4] .+ mixed_count2[:,5]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[:,1],mixed_count3[:,5] ./(mixed_count3[:,4] .+ mixed_count3[:,5]),color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_count4[:,1],mixed_count4[:,5] ./(mixed_count4[:,4] .+ mixed_count4[:,5]),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5)
# plot!(mixed_count5[:,1],mixed_count5[:,5] ./(mixed_count5[:,4] .+ mixed_count5[:,5]),color = :orange3,linestyle =:dash,linewidth=2,
# labels=false,grid=false)
# scatter!(mixed_count5[1:5:end,1],mixed_count5[1:5:end,5] ./(mixed_count5[1:5:end,4] .+ mixed_count5[1:5:end,5]),color = :orange3,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)
# plot!(mixed_count6[1:end,1],mixed_count6[1:end,5] ./(mixed_count6[1:end,4] .+ mixed_count6[1:end,5]),color = :purple,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5)
# scatter!(mixed_count6[3:5:end,1],mixed_count6[3:5:end,5] ./(mixed_count6[3:5:end,4] .+ mixed_count6[3:5:end,5]),color = :purple,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :square,markersize = :1.5)
# yticks=0:0.1:0.4,ylims=[-0.02,0.4]

## Ice mass
# plot(mixed_count1[1:end,1],mixed_count1[1:end,5] .*1e8,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[1:end,1],mixed_count2[1:end,6] .*1e8,color = :blue4,linestyle =:dash,linewidth=2,
# labels="M-1",grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[1:end,1],mixed_count3[1:end,6] .*1e8,color = :green4,linestyle =:dot,linewidth=2.5,
# labels="H-1",grid=false,thickness_scaling=1.5,xticks=0:8:64,xlims=[0,64])
# plot!(mixed_count4[1:end,1],mixed_count4[1:end,6] .*1e8,color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="L-2",grid=false,thickness_scaling=1.5,yticks=0:1:6,ylims=[-0.1,6])

## Total mass
# plot(mixed_count1[:,1],log.(10,(mixed_count1[:,4] .+ mixed_count1[:,5])),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[:,1],log.(10,(mixed_count2[:,4] .+ mixed_count2[:,5])),color = :blue4,linestyle =:dash,linewidth=2,
# labels="H-H",grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[:,1],log.(10,(mixed_count3[:,4] .+ mixed_count3[:,5])),color = :green4,linestyle =:dot,linewidth=2.5,
# labels="L-L",grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_count4[:,1],log.(10,(mixed_count4[:,4] .+ mixed_count4[:,5])),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="H-L",grid=false,thickness_scaling=1.5)
# plot!(mixed_count5[:,1],log.(10,(mixed_count5[:,4] .+ mixed_count5[:,5])),color = :orange3,linestyle =:dash,linewidth=2,
# labels="M-2",grid=false,thickness_scaling=1.5,legend=false)
# scatter!(mixed_count5[1:5:end,1],log.(10,(mixed_count5[1:5:end,4] .+ mixed_count5[1:5:end,5])),color = :orange3,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)

## Total area
# plot(mixed_count1[:,1],log.(10,(mixed_count1[:,6] .+ mixed_count1[:,7])),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[:,1],log.(10,(mixed_count2[:,6] .+ mixed_count2[:,7])),color = :blue4,linestyle =:dash,linewidth=2,
# labels="H-H",grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[:,1],log.(10,(mixed_count3[:,6] .+ mixed_count3[:,7])),color = :green4,linestyle =:dot,linewidth=2.5,
# labels="L-L",grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_count4[:,1],log.(10,(mixed_count4[:,6] .+ mixed_count4[:,7])),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="H-L",grid=false,thickness_scaling=1.5)
# plot!(mixed_count5[:,1],log.(10,(mixed_count5[:,6] .+ mixed_count5[:,7])),color = :orange3,linestyle =:dash,linewidth=2,
# labels="M-2",grid=false,thickness_scaling=1.5,legend=false)
# scatter!(mixed_count5[1:5:end,1],log.(10,(mixed_count5[1:5:end,6] .+ mixed_count5[1:5:end,7])),color = :orange3,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)


## Liquid area
# plot(mixed_count1[1:end,1],mixed_count1[1:end,6] .*1e5,color = :red4,linestyle =:solid,linewidth=2,
# labels="L-1",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

## Ice area
# plot(mixed_count1[1:end,1],mixed_count1[1:end,7] .*1e5,color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


## Relative dipsersions
# b1 = sqrt.(vapor1[:,3])./abs.(vapor1[:,2])
# b2 = sqrt.(vapor2[:,3])./abs.(vapor2[:,2])
# b3 = sqrt.(vapor3[:,3])./abs.(vapor3[:,2])
# b4 = sqrt.(vapor4[:,3])./abs.(vapor4[:,2])
# b5 = sqrt.(vapor5[:,3])./abs.(vapor5[:,2])
# b6 = sqrt.(vapor6[:,3])./abs.(vapor6[:,2])

# b1 = mixed_rad1[1:end,5]
# b2 = mixed_rad2[1:end,5]
# b3 = mixed_rad3[1:end,5]
# b4 = mixed_rad4[1:end,5]
# b5 = mixed_rad5[1:end,5]

# b1 = mixed_rad1[1:end,5] ./ mixed_rad1[1:end,4]
# b2 = mixed_rad2[1:end,5] ./ mixed_rad2[1:end,4]
# b3 = mixed_rad3[1:end,5] ./ mixed_rad3[1:end,4]
# b4 = mixed_rad4[1:end,5] ./ mixed_rad4[1:end,4]
# b5 = mixed_rad5[1:end,5] ./ mixed_rad5[1:end,4]


# # yticks=-1.8:0.4:-0.5,ylims=[-1.8,-0.5]
# plot(mixed_count1[1:end,1],log.(10,b1[1:end]),color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[1:end,1],log.(10,b2[1:end]),color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[1:end,1],log.(10,b3[1:end]),color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_count4[1:end,1],log.(10,b4[1:end]),color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(mixed_count5[1:end,1],log.(10,b5[1:end]),color = :orange3,linestyle =:dash,linewidth=2,
# labels=false)
# scatter!(mixed_count5[1:5:end,1],log.(10,b5[1:5:end]),color = :orange3,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)

# plot(mixed_count1[1:end,1],b1[1:end],color = :red4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(mixed_count2[1:end,1],b2[1:end],color = :blue4,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,legend=false)
# plot!(mixed_count3[1:end,1],b3[1:end],color = :green4,linestyle =:dot,linewidth=2.5,
# labels=false,grid=false,thickness_scaling=1.5,xticks=0:20:120,xlims=[0,120])
# plot!(mixed_count4[1:end,1],b4[1:end],color = :yellow4,linestyle =:dashdot,linewidth=2,
# labels="false",grid=false,thickness_scaling=1.5)
# plot!(mixed_count5[1:end,1],b5[1:end],color = :orange3,linestyle =:dash,linewidth=2,labels=false)
# scatter!(mixed_count5[1:5:end,1],b5[1:5:end],color = :orange3,
# labels=false,grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)

# # Kolmogorov scale
# plot(transition1[:,1],transition1[:,3],color = :red4,linestyle =:solid,linewidth=2,
# labels="L-H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(transition2[:,1],transition2[:,3],color = :black,linestyle =:dashdot,linewidth=2,
# labels="H-H",grid=false,thickness_scaling=1.5,markershape = :circle,markersize = :2)
# plot!(transition3[:,1],transition3[:,3],color = :bwr,linestyle =:dot,linewidth=2.5,
# labels="L-L",grid=false,thickness_scaling=1.5,markershape = :square,markersize = :1.5)
# plot!(transition4[:,1],transition4[:,3],color = :black,linestyle =:dashdot,linewidth=2,
# labels="H-L",grid=false,thickness_scaling=1.5)


# # Stokes number
# L = 0.512
# rho_w = 1000
# rho_a = 1
# mu = 1.5e-5
# g = 9.8

# r = 0.00001076238110
# # eps = 0.00557797206054
# eps = 0.00060994446801
# tau_p = (2 *rho_w*(r^2))/(9*mu)
# nu = mu/rho_a
# tau_f = (nu/eps)^0.5
# St = tau_p/tau_f


# Legends in Plots
# r = zeros(10)
# R = zeros(10)

# plot(r,R,color = :red4,linestyle =:dot,linewidth=2.5,
# labels="L, σ=0.005",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :blue4,linestyle =:dash,linewidth=2,
# labels="L-L, σ=0.13",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :green4,linestyle =:dashdot,linewidth=2,
# labels="L-H, σ=0.11",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :orange3,linestyle =:solid,linewidth=2,
# labels="H, σ=0.0008",grid=false,thickness_scaling=1.5,ickfont=font(9,"Helvetica Bold"))
# plot!(r,R,color = :black,linestyle =:dash,linewidth=2,markershape = :square,markersize = :2,
# labels="H-L, σ=0.10",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :purple,linestyle =:dashdot,linewidth=2,markershape = :circle,markersize = :1.5,
# labels="L-LR",grid=false,legendfont=font(5,"Helvetica Bold"))

# plot(r,R,color = :summer,linestyle =:solid,linewidth=2,
# labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :bwr,linestyle =:dash,linewidth=2,
# labels="L",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :purple,linestyle =:dot,linewidth=2.5,
# labels="H0",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))

# plot(r,R,color = :purple,linestyle =:dot,linewidth=2.5,
# labels="τ_L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(r,R,color = :blue4,linestyle =:solid,linewidth=2,
# labels="τ_λ",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :orange3,linestyle =:dash,linewidth=2,markershape = :circle,markersize = :2,
# labels="τ_η",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
# plot!(r,R,color = :green4,linestyle =:dash,linewidth=2,legendfont=font(5,"Helvetica Bold"),
# labels="τ_m",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(r,R,color = :red4,linestyle =:solid,linewidth=2,markershape = :square,markersize = :2,
# labels="τ_B",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(r,R,color = :black,linestyle =:dashdot,linewidth=2,
# labels="τ_p",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# plot(r,R,color = :purple,linestyle =:dash,linewidth=2,legendfont=font(5,"Helvetica Bold"),
# labels="Case L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

#=
P1 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-radius/radius-0.00", skipstart=3)
P2 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-radius/radius-2.00", skipstart=2)
P3 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-radius/radius-6.00", skipstart=2)

plot(P1[:,1].*1e6,log.(10,P1[:,2]),color = :red4,linestyle =:dot,linewidth=2.5,
labels="false",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(P2[:,1].*1e6,log.(10,P2[:,2]),color = :blue4,linestyle =:dash,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5,xticks=0:0.5:3,xlims=[0,3],yticks=0:1:8,ylims=[0.5,8])
plot!(P3[:,1].*1e6,log.(10,P3[:,2]),color = :green4,linestyle =:solid,linewidth=2,
labels="false",grid=false,thickness_scaling=1.5,legend=false)


P1 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-supersat/supersat-0.00", skipstart=3)
P2 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-supersat/supersat-2.00", skipstart=2)
P3 = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/PDF-supersat/supersat-6.00", skipstart=2)

# plot(P1[:,1],P1[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="0 s",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(P2[:,1],P2[:,2],color = :blue4,linestyle =:dash,linewidth=2,
# labels="2 s",grid=false,thickness_scaling=1.5,yticks=0:1000:3000,ylims=[-100,3000])
# plot!(P3[:,1],P3[:,2],color = :green4,linestyle =:solid,linewidth=2,
# labels="6 s",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"),legend=:top)

plot(P1[:,1],log.(10,P1[:,2]),color = :red4,linestyle =:dot,linewidth=2.5,
labels="0 s",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
plot!(P2[:,1],log.(10,P2[:,2]),color = :blue4,linestyle =:dash,linewidth=2,
labels="2 s",grid=false,thickness_scaling=1.5)
plot!(P3[:,1],log.(10,P3[:,2]),color = :green4,linestyle =:solid,linewidth=2,
labels="6 s",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"),legend=:top)
=#


#=
# One-point correlation for PoF/cloud
RN1 = readdlm("PR_DNS_Scalar/F3_CaseH2_256/RN", skipstart=1)
supersat1 = readdlm("PR_DNS_Scalar/F3_CaseH2_256/supersat")
RN2 = readdlm("PR_DNS_Scalar/F3_CaseL2_256/RN", skipstart=1)
supersat2 = readdlm("PR_DNS_Scalar/F3_CaseL2_256/supersat")

a1 = RN1[1:49,1]
b1 = RN1[1:49,4].*sqrt.(supersat1[1:49,3]).*1e6
c1 = RN1[1:49,3].*supersat1[1:49,2].*1e6
a2 = RN2[1:49,1]
b2 = RN2[1:49,4].*sqrt.(supersat2[1:49,3]).*1e6
c2 = RN2[1:49,3].*supersat2[1:49,2].*1e6

plot(a1,b1./c1,color = :yellow4,linestyle =:solid,linewidth=2,
labels="H-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

plot!(a2,b2./c2,color = :blue4,linestyle =:dash,linewidth=2,
labels="M-L",grid=false,thickness_scaling=1.5,legend=false)

# # Validation
# A1 = readdlm("Post_Process_Scalar/Data/G_2018.txt",skipstart=1)
# A2 = readdlm("Post_Process_Scalar/Data/G_256.txt",skipstart=1)

# plot(A1[:,1],A1[:,2],color = :purple,linestyle =:dash,linewidth=2,yticks=-50:10:-10,ylims=[-50,-10],
# labels="H-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot!(A2[:,1],A2[:,2],color = :bwr,linestyle =:dot,linewidth=2.5,
# labels="M-L",grid=false,thickness_scaling=1.5,legend=false)

# Activation fraction
A = readdlm("Post_Process_Scalar/Data/activ_H.txt",skipstart=1)
A[1,1] = 0.0
A[2:end,1] .= A[2:end,1] .+ 0.2
A[1:4,2] .= 0.0
plot(A[:,1],A[:,2],color = :yellow4,linestyle =:solid,linewidth=2,legend=false,
labels="V",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

A = readdlm("Post_Process_Scalar/Data/activ_H_M.txt",skipstart=1)
A[1,1] = 0.0
A[2:end,1] .= A[2:end,1] .+ 0.2
A[1:4,2] .= 0.0
plot!(A[:,1],A[:,2],color = :black,linestyle =:dot,linewidth=2.5,legend=false,
labels="V",grid=false,thickness_scaling=1.5,xticks=0:5:25,xlims=[0,25])

A = readdlm("Post_Process_Scalar/Data/activ_H_L.txt",skipstart=1)
A[1,1] = 0.0
A[2:end,1] .= A[2:end,1] .+ 0.2
A[1:4,2] .= 0.0
plot!(A[:,1],A[:,2],color = :green4,linestyle =:dash,linewidth=2,legend=false,
labels="V",grid=false,thickness_scaling=1.5)
=#

# vapor = readdlm("PR_DNS_Scalar/F3_CaseL1_256_F0/vapor")
# A1 = vapor1[:,1]
# A2 = zeros(81)

# A2[1:21] = vapor[1:21,2]
# A2[22:81] .= vapor[21,2]

# plot(A1,A2,color = :bwr,linestyle =:dash,linewidth=2,yticks=3.94:0.01:3.98,ylims=[3.94,3.98],legend=false,
# labels="H-L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


#=
A = readdlm("Post_Process_Scalar/particle_count/L2_256.txt")
np = 16e6
plot(A[:,1],A[:,2]./np,color = :red4,linestyle =:dot,linewidth=2.5,legend=false,
labels="V",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

A = readdlm("Post_Process_Scalar/particle_count/L2_256_F.txt",skipstart=1)
plot!(A[:,1],A[:,2],color = :blue4,linestyle =:dash,linewidth=2,legend=false,
labels="V",grid=false,xticks=0:5:40)

A = readdlm("Post_Process_Scalar/particle_count/H2_256.txt")
plot!(A[:,1],A[:,2]./np,color = :yellow4,linestyle =:dashdot,linewidth=2,legend=false,
labels="V",grid=false)

A = readdlm("Post_Process_Scalar/particle_count/H2_256_F.txt",skipstart=1)
plot!(A[1:end,1],A[1:end,2],color = :black,linestyle =:solid,linewidth=2,legend=false,
labels="V",grid=false)
=#



# # Upscale Transfer of Energy and Thermal Boundary
# RN1 = readdlm("PR_DNS_Upscale/Case_H_P_256/RN", skipstart=2)
# supersat1 = readdlm("PR_DNS_Upscale/Case_H_P_256/supersat")
# temperature1 = readdlm("PR_DNS_Upscale/Case_H_P_256/temperature")
# transition1 = readdlm("PR_DNS_Upscale/Case_H_P_256/transition", skipstart=2)
# vapor1 = readdlm("PR_DNS_Upscale/Case_H_P_256/vapor")
# Cd1 = readdlm("PR_DNS_Upscale/Case_H_P_256/Cd", skipstart=1)

# RN2 = readdlm("PR_DNS_Upscale/Case_L_P_256/RN", skipstart=2)
# supersat2 = readdlm("PR_DNS_Upscale/Case_L_P_256/supersat")
# temperature2 = readdlm("PR_DNS_Upscale/Case_L_P_256/temperature")
# transition2 = readdlm("PR_DNS_Upscale/Case_L_P_256/transition", skipstart=2)
# vapor2 = readdlm("PR_DNS_Upscale/Case_L_P_256/vapor")
# Cd2 = readdlm("PR_DNS_Upscale/Case_L_P_256/Cd", skipstart=1)


# RN2 = readdlm("PR_DNS_Upscale/Case_H_D_128/RN", skipstart=2)
# supersat2 = readdlm("PR_DNS_Upscale/Case_H_D_128/supersat")
# temperature2 = readdlm("PR_DNS_Upscale/Case_H_D_128/temperature")
# transition2 = readdlm("PR_DNS_Upscale/Case_H_D_128/transition", skipstart=2)
# vapor2 = readdlm("PR_DNS_Upscale/Case_H_D_128/vapor")
# Cd2 = readdlm("PR_DNS_Upscale/Case_H_D_128/Cd", skipstart=1)

# plot(temperature1[:,1],temperature1[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(temperature2[:,1],temperature2[:,2],color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(vapor1[:,1],vapor1[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(vapor2[:,1],vapor2[:,2],color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(supersat1[:,1],supersat1[:,2],color = :summer,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],supersat2[:,2],color = :bwr,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20))

# plot(supersat1[:,1],sqrt.(supersat1[:,3]),color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],sqrt.(supersat2[:,3]),color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20))

# plot(supersat1[:,1],sqrt.(supersat1[:,3]) ./ supersat1[:,2],color = :summer,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],sqrt.(supersat2[:,3]) ./ supersat2[:,2],color = :bwr,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20),yticks=-25:5:0,ylims=(-26,0))

# plot(RN1[1:37,1],RN1[1:37,3] .*1e6,color = :summer,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[1:37,1],RN2[1:37,3] .*1e6,color = :bwr,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20))

# plot(RN1[1:37,1],RN1[1:37,4] .*1e6,color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[1:37,1],RN2[1:37,4] .*1e6,color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20))

# plot(RN1[1:37,1],RN1[1:37,4] ./RN1[1:37,3],color = :summer,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(RN2[1:37,1],RN2[1:37,4] ./RN2[1:37,3],color = :bwr,linestyle =:dash,linewidth=2,
# labels=false,grid=false,xticks=0:5:20,xlims=(0,20))


## Relative dipsersions
# b1 = sqrt.(temperature1[:,3])./abs.(temperature1[:,2])
# b2 = sqrt.(temperature2[:,3])./abs.(temperature2[:,2])

# b1 = sqrt.(vapor1[:,3])./abs.(vapor1[:,2])
# b2 = sqrt.(vapor2[:,3])./abs.(vapor2[:,2])

# b1 = sqrt.(supersat1[:,3])./abs.(supersat1[:,2])
# b2 = sqrt.(supersat2[:,3])./abs.(supersat2[:,2])

# b1 = RN1[1:end,4] ./ RN1[1:end,3]
# b2 = RN2[1:end,4] ./ RN2[1:end,3]

# plot(supersat1[:,1],log.(10,b1),color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(supersat2[:,1],log.(10,b2),color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

## TKE and activation/deactivation statistics
# plot(transition1[:,1],transition1[:,2],color = :yellow4,linestyle =:solid,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# plot!(transition2[:,1],transition2[:,2],color = :nuuk,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


