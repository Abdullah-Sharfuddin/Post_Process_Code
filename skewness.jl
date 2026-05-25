using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

function nodal_value(var,scalar,Nx,Ny,Nz)

    n = 192
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
N = 128
Nx = 8
Ny = 4
Nz = 4

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,4])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuation
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Normalize
mean = sum(scalar[:])/((n-2)^3)
fluctuation = scalar[:] .- mean
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms


# Compute skewness of normalized scalar
avg = sum(scalar[:])/((n-2)^3)
prime = scalar[:] .- avg

first_term = sum(prime[:].^3) / (n-2)^3
second_term = sum(prime[:].^2) / (n-2)^3
second_term = second_term^(3/2)

skewness = first_term/second_term

