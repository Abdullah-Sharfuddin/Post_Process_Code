using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

function nodal_value(var,scalar,Nx,Ny,Nz)

    n = 288
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


var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_P/nodal_values/nodal_values-4.00-$I")
    append!(var,VAR[:,9])
end

n = 288
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# sc = scalar[:]
# a = 0
# for m = 1:n^3
#     global a
#     if abs(sc[m]) > 10.0
#         a = a + 1
#     end
# end


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
fluctuation = scalar[:] .- mean_scalar
rms = sum(fluctuation[:].^2)/((n-2)^3)
rms = sqrt(rms)
scalar = scalar ./rms

# Compute kurtosis
avg = sum(scalar[:])/((n-2)^3)
prime = scalar[:] .- avg

first_term = sum(prime[:].^4) / (n-2)^3
second_term = sum(prime[:].^2) / (n-2)^3
second_term = second_term^(4/2)

kurtosis = first_term/second_term
