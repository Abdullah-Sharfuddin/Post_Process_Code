
using DelimitedFiles
using Distributions
using PlotlyJS
# using Plots
using ColorSchemes
using LaTeXStrings
using LinearAlgebra


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

    return scalar
end


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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_L/nodal_values/nodal_values-4.00-$I")
    append!(var,VAR[:,4])
end

n = 128
scalar = zeros(n,n,n)
scalar = nodal_value(var,scalar,Nx,Ny,Nz)

# # Read data
# N = 128
# Nx = 8
# Ny = 4
# Nz = 4

# var = []
# I = 0
# for i = 1:N
#     global I, var
#     I = i - 1
#     VAR = readdlm("PR_DNS_Upscale/Case_H_P_128/nodal_values/nodal_values-10.00-$I")
#     append!(var,VAR[:,3])
# end

# n = 128
# scalar = zeros(n,n,n)
# scalar = nodal_value(var,scalar,Nx,Ny,Nz)

# # Compute mean field
# mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuating field
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#             scalar[i,j,k] = scalar[i,j,k]/mean_scalar
#         end
#     end
# end

L = 0.512
h = L/(n-1)

# X = collect(0:h:L)
# plot(X,scalar[:,64,64], color = :purple,linestyle =:solid,linewidth=2,tickfont=font(9,"Helvetica Bold"),
# labels="H-L",grid=false,thickness_scaling=1.5,legend=false)

# # Normalize
# mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
# std = 0.0
# for k=1:n
#     for j=1:n
#         for i=1:n
#         global std
#         std = std + sqrt((scalar[i,j,k]-mean_scalar)^2)
#         end
#     end
# end

# std = std / n^3
# scalar[:] .= scalar[:] .- mean_scalar
# scalar[:] .= scalar[:] ./ std

# # Plot for 128^3 points # The PlotlyJS library cannot handle 256^3 data
# nn = 128
# new_scalar = zeros(nn,nn,nn)
# for k=1:nn
#     for j=1:nn
#         for i=1:nn
#             m = 2*i - 1
#             new_scalar[i,j,k] = scalar[m,m,m] - mean_scalar
#         end
#     end
# end


nn = n
data = range(0, stop = L, length=nn)
X, Y, Z = mgrid(data, data, data)

layout = Layout(
    scene=attr(
        xaxis=attr(
            nticks=6,tickfont=attr(size=14,family="Arial",color=:black,),
        ),
        yaxis=attr(
            nticks=6,tickfont=attr(size=14,family="Arial",color=:black,),
        ),
        zaxis=attr(
            nticks=6,tickfont=attr(size=14,family="Arial",color=:black,),
        ),
    ),
)

trace = isosurface(
    x=X[:],y=Y[:],z=Z[:],value=scalar[:],
    colorscale=colors.Spectral,
)

plot(trace,layout)

# plot(isosurface(x=Z[:],y=Y[:],z=X[:],value=new_scalar[:],colorscale=colors.Spectral),layout,)


#=
# Plot the colorbar
using DelimitedFiles
using Plots
using LaTeXStrings

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

    return scalar
end

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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_L/nodal_values/nodal_values-30.00-$I")
    append!(var,VAR[:,6])
end

n = 128
scalar = zeros(n,n,n)
scalar = nodal_value(var,scalar,Nx,Ny,Nz)


pyplot()
contourf(scalar[:,64,:],
    levels=5,
    colorbar_tickfontsize=12,
    # colorbar_ticks=270.750:0.003:270.763,
    # clims=(270.750,270.763),
    c=:Spectral
)
=#
#=
using Plots
# Read data
N = 128
r = []
x = []
y = []
z = []
rho = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rc, rr
    K = i - 1
    if (filesize("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-0.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-0.00-$K")
        append!(r,rr[:,1])
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
r_w = zeros(count_w)
r_ice = zeros(count_ice)
count_1 = 0
count_2 = 0
for i = 1:nn
    global X_w, X_ice, count_1, count_2, r_w, r_ice
	if rho[i] == 1000.0
        count_1 = count_1 + 1
		X_w[count_1] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_w[count_1] = r[i]
	else
        count_2 = count_2 + 1
		X_ice[count_2] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_ice[count_2] = r[i]
	end
end

nc = 256
D = sqrt(0.512^2 + 0.512^2 + 0.512^2)
h = D/(nc-1)
X1 = collect(0:h:D)
X2 = collect(0:h:D)
r1 = zeros(nc)
r2 = zeros(nc)

for j=1:nc-1
    count_3 = 0
    K1 = 0
    for i = 1:count_w
        global r1,X1
        if ((X_w[i] > X1[j]) && (X_w[i] <= X1[j+1]))
            K1 = K1 + r_w[i]
            count_3 = count_3 + 1
        end
        r1[j] = K1/count_3
    end
end

for j=1:nc-1
    count_4 = 0
    K2 = 0
    for i = 1:count_ice
        global r2,X2
        if ((X_ice[i] > X2[j]) && (X_ice[i] <= X2[j+1]))
            K2 = K2 + r_ice[i]
            count_4 = count_4 + 1
        end
        r2[j] = K2/count_4
    end
end

plot(X2[3:end-2],r2[3:end-2] .*1e6,color = :bwr,linestyle =:solid,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))

# plot!(X2[3:end-2],r2[3:end-2] .*1e6,color = :purple,linestyle =:dash,linewidth=2,
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))


# Read data
N = 128
r = []
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
        append!(r,rr[:,1])
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
r_w = zeros(count_w)
r_ice = zeros(count_ice)
count_1 = 0
count_2 = 0
for i = 1:nn
    global X_w, X_ice, count_1, count_2, r_w, r_ice
	if rho[i] == 1000.0
        count_1 = count_1 + 1
		X_w[count_1] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_w[count_1] = r[i]
	else
        count_2 = count_2 + 1
		X_ice[count_2] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_ice[count_2] = r[i]
	end
end

nc = 256
D = sqrt(0.512^2 + 0.512^2 + 0.512^2)
h = D/(nc-1)
X1 = collect(0:h:D)
X2 = collect(0:h:D)
r1 = zeros(nc)
r2 = zeros(nc)

for j=1:nc-1
    count_3 = 0
    K1 = 0
    for i = 1:count_w
        global r1,X1
        if ((X_w[i] > X1[j]) && (X_w[i] <= X1[j+1]))
            K1 = K1 + r_w[i]
            count_3 = count_3 + 1
        end
        r1[j] = K1/count_3
    end
end

for j=1:nc-1
    count_4 = 0
    K2 = 0
    for i = 1:count_ice
        global r2,X2
        if ((X_ice[i] > X2[j]) && (X_ice[i] <= X2[j+1]))
            K2 = K2 + r_ice[i]
            count_4 = count_4 + 1
        end
        r2[j] = K2/count_4
    end
end

plot!(X2[3:end-2],r2[3:end-2] .*1e6,color = :summer,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))


# Read data
N = 128
r = []
x = []
y = []
z = []
rho = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rc, rr
    K = i - 1
    if (filesize("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-60.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-60.00-$K")
        append!(r,rr[:,1])
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
r_w = zeros(count_w)
r_ice = zeros(count_ice)
count_1 = 0
count_2 = 0
for i = 1:nn
    global X_w, X_ice, count_1, count_2, r_w, r_ice
	if rho[i] == 1000.0
        count_1 = count_1 + 1
		X_w[count_1] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_w[count_1] = r[i]
	else
        count_2 = count_2 + 1
		X_ice[count_2] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_ice[count_2] = r[i]
	end
end

nc = 256
D = sqrt(0.512^2 + 0.512^2 + 0.512^2)
h = D/(nc-1)
X1 = collect(0:h:D)
X2 = collect(0:h:D)
r1 = zeros(nc)
r2 = zeros(nc)

for j=1:nc-1
    count_3 = 0
    K1 = 0
    for i = 1:count_w
        global r1,X1
        if ((X_w[i] > X1[j]) && (X_w[i] <= X1[j+1]))
            K1 = K1 + r_w[i]
            count_3 = count_3 + 1
        end
        r1[j] = K1/count_3
    end
end

for j=1:nc-1
    count_4 = 0
    K2 = 0
    for i = 1:count_ice
        global r2,X2
        if ((X_ice[i] > X2[j]) && (X_ice[i] <= X2[j+1]))
            K2 = K2 + r_ice[i]
            count_4 = count_4 + 1
        end
        r2[j] = K2/count_4
    end
end

plot!(X2[3:end-2],r2[3:end-2] .*1e6,color = :purple,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
=#

#=
# Read data
N = 128
r = []
x = []
y = []
z = []
rho = []
N1 = 1
N2 = N
for i = N1:N2
    global K, rc, rr
    K = i - 1
    if (filesize("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-120.00-$K") == 0)
        continue
    else
        rr = readdlm("PR_DNS_Mixed/Case_H_S_L/particle/particle_radius-120.00-$K")
        append!(r,rr[:,1])
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
r_w = zeros(count_w)
r_ice = zeros(count_ice)
count_1 = 0
count_2 = 0
for i = 1:nn
    global X_w, X_ice, count_1, count_2, r_w, r_ice
	if rho[i] == 1000.0
        count_1 = count_1 + 1
		X_w[count_1] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_w[count_1] = r[i]
	else
        count_2 = count_2 + 1
		X_ice[count_2] = sqrt(x[i]^2 +y[i]^2 +z[i]^2)
        r_ice[count_2] = r[i]
	end
end

nc = 256
D = sqrt(0.512^2 + 0.512^2 + 0.512^2)
h = D/(nc-1)
X1 = collect(0:h:D)
X2 = collect(0:h:D)
r1 = zeros(nc)
r2 = zeros(nc)

for j=1:nc-1
    count_3 = 0
    K1 = 0
    for i = 1:count_w
        global r1,X1
        if ((X_w[i] > X1[j]) && (X_w[i] <= X1[j+1]))
            K1 = K1 + r_w[i]
            count_3 = count_3 + 1
        end
        r1[j] = K1/count_3
    end
end

for j=1:nc-1
    count_4 = 0
    K2 = 0
    for i = 1:count_ice
        global r2,X2
        if ((X_ice[i] > X2[j]) && (X_ice[i] <= X2[j+1]))
            K2 = K2 + r_ice[i]
            count_4 = count_4 + 1
        end
        r2[j] = K2/count_4
    end
end

plot!(X2[3:end-2],r2[3:end-2] .*1e6,color = :black,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
=#
