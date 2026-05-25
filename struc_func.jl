using DelimitedFiles
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

# Function to compute the spatial correlation coefficient
function struc_func(v1,v2,v3,v4,n,I)

    # Initialization
    term = zeros(n,n,n)

    # Compute the terms at nodal points
    for k = 1:n
        index_k = k + I
        if index_k > n
            index_k = index_k - n
        end
        for j = 1:n
            index_j = j + I
            if index_j > n
                index_j = index_j - n
            end
            for i = 1:n
                index_i = i + I
                if index_i > n
                    index_i = index_i - n
                end
                term[i,j,k] = (v4[index_i,index_j,index_k]-v4[i,j,k])^3

                # term[i,j,k] = (v4[index_i,index_j,index_k]-v4[i,j,k])^2
                # term[i,j,k] = term[i,j,k] * (v3[index_i,index_j,index_k]-v3[i,j,k])
            end
        end
    end

    # Compute the structure function
    SF = sum(term[:,:,:]) / (n*n*n)

    return SF
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

    return scalar, count
end


# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
var4 = []
I = 0
for i = 1:N
    global I, var1, var2
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,8])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)
scalar3 = zeros(n,n,n)
scalar4 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar1,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var2,scalar2,Nx,Ny,Nz)
(scalar3,count) = nodal_value(var3,scalar3,Nx,Ny,Nz)
(scalar4,count) = nodal_value(var4,scalar4,Nx,Ny,Nz)

# # For liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar4[i,j,k] = C*(scalar4[i,j,k]^3)
#         end
#     end
# end

# # Compute mean field
# mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
# mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)
# mean_scalar3 = sum(scalar3[:,:,:]) / (n*n*n)
# mean_scalar4 = sum(scalar4[:,:,:]) / (n*n*n)

# # Compute fluctuating field
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
#             scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
#             scalar3[i,j,k] = scalar3[i,j,k] - mean_scalar3
#             scalar4[i,j,k] = scalar4[i,j,k] - mean_scalar4
#         end
#     end
# end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
S = zeros(m)

# Calculate the inertial subrange
eta = 0.00088196095823
l = -0.6
u = -1.0
l_limit = (2*pi)/(10^(l))*eta
u_limit = (2*pi)/(10^(u))*eta

# Final loop
for i = 1:m
    ii = i-1
    global S, r

    # Get structure function
    S[i] =  struc_func(scalar1,scalar2,scalar3,scalar4,n,ii)

    # Update separation distance
    r[i] = ii*h

    # Final value
    eps = 0.00557797206054
    S[i] = -S[i]/(eps*r[i])
end

plot(r[1:32],S[1:32],color = :summer,linestyle =:solid,linewidth=2,yticks=0.0005:0.0005:0.0025,ylims=[0.0005,0.0028],
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot(r[1:32],S[1:32] .*1e9,color = :summer,linestyle =:solid,linewidth=2,yticks=-6:1:1,ylims=[-6,1],
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))



# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
var4 = []
I = 0
for i = 1:N
    global I, var1, var2
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,8])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)
scalar3 = zeros(n,n,n)
scalar4 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar1,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var2,scalar2,Nx,Ny,Nz)
(scalar3,count) = nodal_value(var3,scalar3,Nx,Ny,Nz)
(scalar4,count) = nodal_value(var4,scalar4,Nx,Ny,Nz)

# # For liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar4[i,j,k] = C*(scalar4[i,j,k]^3)
#         end
#     end
# end

# # Compute mean field
# mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
# mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)
# mean_scalar3 = sum(scalar3[:,:,:]) / (n*n*n)
# mean_scalar4 = sum(scalar4[:,:,:]) / (n*n*n)

# # Compute fluctuating field
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
#             scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
#             scalar3[i,j,k] = scalar3[i,j,k] - mean_scalar3
#             scalar4[i,j,k] = scalar4[i,j,k] - mean_scalar4
#         end
#     end
# end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
S = zeros(m)

# Calculate the inertial subrange
eta = 0.00153371925558
l = -0.6
u = -1.0
l_limit = (2*pi)/(10^(l))*eta
u_limit = (2*pi)/(10^(u))*eta

# Final loop
for i = 1:m
    ii = i-1
    global S, r

    # Get structure function
    S[i] =  struc_func(scalar1,scalar2,scalar3,scalar4,n,ii)

    # Update separation distance
    r[i] = ii*h

    # Final value
    eps = 0.00060994446801
    S[i] = -S[i]/(eps*r[i])
end

plot!(twinx(),r[1:32],S[1:32] ,color = :bwr,linestyle =:dash,linewidth=2,yticks=-0.01:0.01:0.03,ylims=[-0.015,0.03],
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

# plot!(twinx(),r[1:32],S[1:32]  .*1e8,color = :bwr,linestyle =:dash,linewidth=2,yticks=-10:4:6,ylims=[-10,7],
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


#=
# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
var4 = []
I = 0
for i = 1:N
    global I, var1, var2
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P0_256/nodal_values/nodal_values-8.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,9])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)
scalar3 = zeros(n,n,n)
scalar4 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar1,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var2,scalar2,Nx,Ny,Nz)
(scalar3,count) = nodal_value(var3,scalar3,Nx,Ny,Nz)
(scalar4,count) = nodal_value(var4,scalar4,Nx,Ny,Nz)

# For liquid water mixing ratio
rho_l = 1000
h = 0.512/256
rho_a = 1.0
C = (4*pi*rho_l)/(3*rho_a*(h^3))
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar4[i,j,k] = C*(scalar4[i,j,k]^3)
        end
    end
end

# Compute mean field
mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)
mean_scalar3 = sum(scalar3[:,:,:]) / (n*n*n)
mean_scalar4 = sum(scalar4[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
            scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
            scalar3[i,j,k] = scalar3[i,j,k] - mean_scalar3
            scalar4[i,j,k] = scalar4[i,j,k] - mean_scalar4
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
S = zeros(m)

# Calculate the inertial subrange
eta = 0.00088567694246
l = -0.6
u = -1.0
l_limit = (2*pi)/(10^(l))*eta
u_limit = (2*pi)/(10^(u))*eta

# Final loop
for i = 1:m
    ii = i-1
    global S, r

    # Get structure function
    S[i] =  struc_func(scalar1,scalar2,scalar3,scalar4,n,ii)

    # Update separation distance
    r[i] = ii*h

    # Final value
    eps = 0.00548494684568
    S[i] = -S[i]/(eps*r[i])
end

plot!(r[1:32],S[1:32],color = :purple,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
=#