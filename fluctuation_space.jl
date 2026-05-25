using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

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
#=
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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_L/nodal_values/nodal_values-0.00-$I")
    append!(var,VAR[:,2])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuation
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#         end
#     end
# end


L = 0.512
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

X = collect(0:c:D)

scalar_diag = zeros(n)
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global scalar_diag, count
            if (i == j && j == k)
                count = count + 1
                scalar_diag[count] = scalar[i,j,k]
            end
        end
    end
end

plot(X[1:end],scalar_diag[1:end],color = :bwr,linestyle =:solid,linewidth=2,tickfont=font(9,"Helvetica Bold"),
grid=false,thickness_scaling=1.5,legend=false)

hline!([mean_scalar],color = :bwr,linestyle =:dash,linewidth=2,tickfont=font(9,"Helvetica Bold"),
grid=false)
=#

#=
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
    append!(var,VAR[:,2])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# # Compute mean
# mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuation
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#         end
#     end
# end


L = 0.512
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

X = collect(0:c:D)

scalar_diag = zeros(n)
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global scalar_diag, count
            if (i == j && j == k)
                count = count + 1
                scalar_diag[count] = scalar[i,j,k]
            end
        end
    end
end

plot(X[1:end],scalar_diag[1:end],color = :summer,linestyle =:dash,linewidth=2,tickfont=font(9,"Helvetica Bold"),
grid=false,thickness_scaling=1.5,legend=false)


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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_L/nodal_values/nodal_values-60.00-$I")
    append!(var,VAR[:,2])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# # Compute mean
# mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuation
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#         end
#     end
# end

L = 0.512
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

X = collect(0:c:D)

scalar_diag = zeros(n)
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global scalar_diag, count
            if (i == j && j == k)
                count = count + 1
                scalar_diag[count] = scalar[i,j,k]
            end
        end
    end
end

plot!(X[1:end],scalar_diag[1:end],color = :purple,linestyle =:dashdot,linewidth=2,tickfont=font(9,"Helvetica Bold"),
grid=false,thickness_scaling=1.5,legend=false)


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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_L/nodal_values/nodal_values-120.00-$I")
    append!(var,VAR[:,2])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# # Compute mean
# mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuation
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#         end
#     end
# end


L = 0.512
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

X = collect(0:c:D)

scalar_diag = zeros(n)
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global scalar_diag, count
            if (i == j && j == k)
                count = count + 1
                scalar_diag[count] = scalar[i,j,k]
            end
        end
    end
end

plot!(X[1:end],scalar_diag[1:end],color = :black,linestyle =:dot,linewidth=2.5,tickfont=font(9,"Helvetica Bold"),
grid=false,thickness_scaling=1.5,legend=false)
=#

#=
# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_P_128/nodal_values/nodal_values-30.00-$I")
    append!(var,VAR[:,9])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# a = 0
# for i = 1:n
#     for j = 1:n
#         for k = 1:n
#             global a
#             if scalar[i,j,k] > 267
#                 a = a + 1
#             end
#         end
#     end
# end

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# # Compute fluctuating field
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = scalar[i,j,k] - mean_scalar
#             # scalar[i,j,k] = scalar[i,j,k]/mean_scalar
#         end
#     end
# end


# L = 0.512
# c = L/(n-1)
# X = collect(0:c:L)

L = 0.512
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)
X = collect(0:c:D)

scalar_diag = zeros(n)
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global scalar_diag, count
            if (i == j && j == k)
                count = count + 1
                scalar_diag[count] = scalar[i,j,k]
            end
        end
    end
end

plot(X[1:end],scalar_diag[1:end],color = :summer,linestyle =:solid,linewidth=2,tickfont=font(9,"Helvetica Bold"),
grid=false,thickness_scaling=1.5,legend=false)
=#

#=
# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var = []
I = 0
for i = 1:N
    global I, var1, var2
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_P_128/particle/particle_radius-10.00-$I")
    append!(var1,VAR[:,2])
    append!(var2,VAR[:,4])
end

np = size(var1,1)
XY = [var1 var2]

sortslices(XY, dims=1)

scatter(XY[1:1000,1],XY[1:1000,2],color = :summer,markershape = :circle,markersize = :2,
grid=false,legend=false)
=#
