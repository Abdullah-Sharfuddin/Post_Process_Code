using DelimitedFiles
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

# Function to compute the spatial correlation coefficient
function corr_coeff(vel,scalar,n,I)

    # Initialization
    term_1 = zeros(n,n,n)
    term_2 = zeros(n,n,n)
    term_3 = zeros(n,n,n)

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
                term_1[i,j,k] = (scalar[i,j,k]*vel[index_i,index_j,index_k])
                term_2[i,j,k] = scalar[i,j,k]^2
                term_3[i,j,k] = vel[i,j,k]^2
            end
        end
    end

    # Compute average
    avg_term_1 = sum(term_1[:,:,:]) / (n*n*n)
    avg_term_2 = sum(term_2[:,:,:]) / (n*n*n)
    avg_term_3 = sum(term_3[:,:,:]) / (n*n*n)

    # Compute the correlation coefficient
    RR = avg_term_1/sqrt((avg_term_2*avg_term_3))

    return RR
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

#=
# Read data
np = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:np
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,5])
    append!(var_2,VAR[:,5])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end
# R .= -R
# R[56:96] .= R[56]

# A = readdlm("Post_Process_Scalar/Data/qv_S_L.txt",skipstart=1)
# plot(A[:,1],A[:,2],color = :red4,linestyle =:dot,linewidth=2.5,
# labels="M",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))

plot(r,R,color = :red4,linestyle =:dot,linewidth=2.5,
labels="M",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# yticks=-0.5:0.25:0.5,ylims=[-0.5,0.5],


# Read data
np = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:np
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,3])
    append!(var_2,VAR[:,3])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end
# R[60:96] .= R[60]
plot!(r,R,color = :blue4,linestyle =:dash,linewidth=2,legend=false,
labels="H",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))



# Read data
np = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:np
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_H_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,3])
    append!(var_2,VAR[:,3])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end
# R[60:96] .= R[60]
plot!(r,R,color = :green4,linestyle =:dashdot,linewidth=2,
labels="V",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,3])
    append!(var_2,VAR[:,3])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end
# R .= -R
# R[56:96] .= R[56]

# A = readdlm("Post_Process_Scalar/Data/qv_S_H.txt",skipstart=1)
# plot!(A[:,1],A[:,2],color = :orange3,linestyle =:solid,linewidth=2,
# labels="V",grid=false,thickness_scaling=1.5)

plot!(r,R,color = :orange3,linestyle =:solid,linewidth=2,
labels="V",grid=false,thickness_scaling=1.5)


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,3])
    append!(var_2,VAR[:,3])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end
# R[60:96] .= R[60]
plot!(r,R,color = :black,linestyle =:dash,linewidth=2,
labels="V",grid=false,thickness_scaling=1.5)
scatter!(r[1:5:end],R[1:5:end],color = :black,markershape = :square,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)
=#

#=
# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var_1 = []
var_2 = []
var_3 = []
I = 0

for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_LA_1/nodal_values/nodal_values-25.00-$I")
    append!(var_1,VAR[:,5])
    append!(var_2,VAR[:,5])
    append!(var_3,VAR[:,8])
end
# var_2 .= var_2 .*(-1)

n = 192
scalar = zeros(n,n,n)
vel = zeros(n,n,n)
forcing = zeros(n,n,n)

(scalar,count) = nodal_value(var_1,scalar,Nx,Ny,Nz)
(vel,count) = nodal_value(var_2,vel,Nx,Ny,Nz)
(forcing,count) = nodal_value(var_3,forcing,Nx,Ny,Nz)
# vel .= vel .+ forcing

# Compute mean field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)
mean_vel = sum(vel[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
            vel[i,j,k] = vel[i,j,k] - mean_vel
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(vel,scalar,n,ii)

    # Update separation distance
    r[i] = ii*h
end

plot!(r,R,color = :purple,linestyle =:dashdot,linewidth=2,
labels="V",grid=false,thickness_scaling=1.5)
scatter!(r[1:5:end],R[1:5:end],color = :purple,markershape = :circle,markersize = :2,
labels=false,grid=false,thickness_scaling=1.5)
=#


## Upscale
# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,9])
    append!(var2,VAR[:,6])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar1,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var2,scalar2,Nx,Ny,Nz)

# For liquid water mixing ratio
rho_l = 1000
h = 0.512/256
rho_a = 1.0
C = (4*pi*rho_l)/(3*rho_a*(h^3))
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar1[i,j,k] = C*(scalar1[i,j,k]^3)
            # scalar2[i,j,k] = C*(scalar2[i,j,k]^3)
        end
    end
end

# Compute mean field
mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
            scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(scalar1,scalar2,n,ii)

    # Update separation distance
    r[i] = ii*h
end

plot(r,R,color = :summer,linestyle =:solid,linewidth=2,yticks=-0.25:0.25:1,ylims=[-0.25,1],
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))



# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,9])
    append!(var2,VAR[:,6])
end

n = 256
scalar1 = zeros(n,n,n)
scalar2 = zeros(n,n,n)

(scalar1,count) = nodal_value(var1,scalar1,Nx,Ny,Nz)
(scalar2,count) = nodal_value(var2,scalar2,Nx,Ny,Nz)

# For liquid water mixing ratio
rho_l = 1000
h = 0.512/256
rho_a = 1.0
C = (4*pi*rho_l)/(3*rho_a*(h^3))
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar1[i,j,k] = C*(scalar1[i,j,k]^3)
            # scalar2[i,j,k] = C*(scalar2[i,j,k]^3)
        end
    end
end

# Compute mean field
mean_scalar1 = sum(scalar1[:,:,:]) / (n*n*n)
mean_scalar2 = sum(scalar2[:,:,:]) / (n*n*n)

# Compute fluctuating field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar1[i,j,k] = scalar1[i,j,k] - mean_scalar1
            scalar2[i,j,k] = scalar2[i,j,k] - mean_scalar2
        end
    end
end

# Spatial separation
L = 0.512
h = L/n
m = Int(floor(n*0.5))
r = zeros(m)

# Coefficient variables
R = zeros(m)

# Final loop
for i = 1:m
    ii = i-1
    global R, r

    # Get correlation coefficient
    R[i] =  corr_coeff(scalar1,scalar2,n,ii)

    # Update separation distance
    r[i] = ii*h
end

plot!(r,R,color = :bwr,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))



