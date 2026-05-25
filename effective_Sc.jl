using DelimitedFiles
using Distributions
using Plots
# using PlotlyJS
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

var1 = []
var2 = []
var3 = []
var4 = []
var5 = []
var6 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,5])
    append!(var5,VAR[:,7])
    append!(var6,VAR[:,8])
end

n = 192
u = zeros(n,n,n)
v = zeros(n,n,n)
w = zeros(n,n,n)
scalar = zeros(n,n,n)
forcing = zeros(n,n,n)
latent = zeros(n,n,n)

(u,count) = nodal_value(var1,u,Nx,Ny,Nz)
(v,count) = nodal_value(var2,v,Nx,Ny,Nz)
(w,count) = nodal_value(var3,w,Nx,Ny,Nz)
(scalar,count) = nodal_value(var4,scalar,Nx,Ny,Nz)
(forcing,count) = nodal_value(var5,forcing,Nx,Ny,Nz)
(latent,count) = nodal_value(var6,latent,Nx,Ny,Nz)

# Calculate the advection term
L = 0.512
h = L/(n-1)

advec = zeros(n,n,n)
for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            global advec
            advec[i,j,k] = u[i,j,k]*((scalar[i+1,j,k] - scalar[i,j,k])/h) + v[i,j,k]*((scalar[i,j+1,k] - scalar[i,j,k])/h) + w[i,j,k]*((scalar[i,j,k+1] - scalar[i,j,k])/h)
        end
    end
end
advec[n,n,n] = u[n,n,n]*((scalar[1,n,n] - scalar[n,n,n])/h) + v[n,n,n]*((scalar[n,1,n] - scalar[n,n,n])/h) + w[n,n,n]*((scalar[n,n,1] - scalar[n,n,n])/h)

# Calculate the Laplacian of the scalar
lap = zeros(n,n,n)
for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global lap
            lap[i-1,j-1,k-1] = (scalar[i+1,j,k] - 2*scalar[i,j,k] + scalar[i-1,j,k])/(h^2) + (scalar[i,j+1,k] - 2*scalar[i,j,k] + scalar[i,j-1,k])/(h^2) + (scalar[i,j,k+1] - 2*scalar[i,j,k] + scalar[i,j,k-1])/(h^2)
        end
    end
end

lap[1,1,1] = (scalar[2,1,1] - 2*scalar[1,1,1] + scalar[n,1,1])/(h^2) + (scalar[1,2,1] - 2*scalar[1,1,1] + scalar[1,n,1])/(h^2) + (scalar[1,1,2] - 2*scalar[1,1,1] + scalar[1,1,n])/(h^2)

lap[n,n,n] = (scalar[1,n,n] - 2*scalar[n,n,n] + scalar[n-1,n,n])/(h^2) + (scalar[n,1,n] - 2*scalar[n,n,n] + scalar[n,n-1,n])/(h^2) + (scalar[n,n,1] - 2*scalar[n,n,n] + scalar[n,n,n-1])/(h^2)


# Calculate artificial forcing viscosity, TKE, and scalar variance
mu_eff = zeros(n,n,n)
TKE = zeros(n,n,n)
mu_v = 2.2e-5
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff, TKE, count
            if (lap[i,j,k] != 0)
                mu_eff[i,j,k] = (advec[i,j,k] - latent[i,j,k] - forcing[i,j,k])/lap[i,j,k]
            end
            TKE[i,j,k] = 0.5*(u[i,j,k]^2 + v[i,j,k]^2 + w[i,j,k]^2)
        end
    end
end
#=
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

mu_eff_diag = zeros(n)
TKE_diag = zeros(n)
forcing_diag = zeros(n)
latent_diag = zeros(n)
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff_diag, count
            if (i == j && j == k)
                count = count + 1
                mu_eff_diag[count] = mu_eff[i,j,k] 
                TKE_diag[count] = TKE[i,j,k]
                forcing_diag[count] = forcing[i,j,k] 
                latent_diag[count] = latent[i,j,k]  
            end
        end
    end
end

X = collect(0:c:D)

plot(X[1:end],mu_eff_diag[1:end],color = :red4,linestyle =:dot,linewidth=2.5,tickfont=font(9,"Helvetica Bold"),
labels="H-L",grid=false,thickness_scaling=1.5,legend=false)
=#

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(mu_eff))/(N^3)
Ek = scalar_hat.^2
Ek = fftshift(Ek)

box_sidex = N
box_sidey = N
box_sidez = N

box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

centerx = Int(box_sidex/2)
centery = Int(box_sidey/2)
centerz = Int(box_sidez/2)

Ek_avsphr = zeros(box_radius-1)
wavenumber = zeros(box_radius-1)

for k = 1:N
	for j = 1:N
		for i = 1:N
            global m
			m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
            if (m <= box_radius-1)
                Ek_avsphr[m] = Ek_avsphr[m] + Ek[i,j,k]
                wavenumber[m] = (2*pi*m/L)
            end
        end
    end
end

plot(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :red4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var3 = []
var4 = []
var5 = []
var6 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,5])
    append!(var5,VAR[:,7])
    append!(var6,VAR[:,8])
end

n = 192
u = zeros(n,n,n)
v = zeros(n,n,n)
w = zeros(n,n,n)
scalar = zeros(n,n,n)
forcing = zeros(n,n,n)
latent = zeros(n,n,n)

(u,count) = nodal_value(var1,u,Nx,Ny,Nz)
(v,count) = nodal_value(var2,v,Nx,Ny,Nz)
(w,count) = nodal_value(var3,w,Nx,Ny,Nz)
(scalar,count) = nodal_value(var4,scalar,Nx,Ny,Nz)
(forcing,count) = nodal_value(var5,forcing,Nx,Ny,Nz)
(latent,count) = nodal_value(var6,latent,Nx,Ny,Nz)


# Calculate the advection term
L = 0.512
h = L/(n-1)

advec = zeros(n,n,n)
for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            global advec
            advec[i,j,k] = u[i,j,k]*((scalar[i+1,j,k] - scalar[i,j,k])/h) + v[i,j,k]*((scalar[i,j+1,k] - scalar[i,j,k])/h) + w[i,j,k]*((scalar[i,j,k+1] - scalar[i,j,k])/h)
        end
    end
end
advec[n,n,n] = u[n,n,n]*((scalar[1,n,n] - scalar[n,n,n])/h) + v[n,n,n]*((scalar[n,1,n] - scalar[n,n,n])/h) + w[n,n,n]*((scalar[n,n,1] - scalar[n,n,n])/h)

# Calculate the Laplacian of the scalar
lap = zeros(n,n,n)
for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global lap
            lap[i-1,j-1,k-1] = (scalar[i+1,j,k] - 2*scalar[i,j,k] + scalar[i-1,j,k])/(h^2) + (scalar[i,j+1,k] - 2*scalar[i,j,k] + scalar[i,j-1,k])/(h^2) + (scalar[i,j,k+1] - 2*scalar[i,j,k] + scalar[i,j,k-1])/(h^2)
        end
    end
end

lap[1,1,1] = (scalar[2,1,1] - 2*scalar[1,1,1] + scalar[n,1,1])/(h^2) + (scalar[1,2,1] - 2*scalar[1,1,1] + scalar[1,n,1])/(h^2) + (scalar[1,1,2] - 2*scalar[1,1,1] + scalar[1,1,n])/(h^2)

lap[n,n,n] = (scalar[1,n,n] - 2*scalar[n,n,n] + scalar[n-1,n,n])/(h^2) + (scalar[n,1,n] - 2*scalar[n,n,n] + scalar[n,n-1,n])/(h^2) + (scalar[n,n,1] - 2*scalar[n,n,n] + scalar[n,n,n-1])/(h^2)


# Calculate artificial forcing viscosity, TKE, and scalar variance
mu_eff = zeros(n,n,n)
TKE = zeros(n,n,n)
mu_v = 2.2e-5
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff, TKE, count
            if (lap[i,j,k] != 0)
                mu_eff[i,j,k] = (advec[i,j,k] - latent[i,j,k] - forcing[i,j,k])/lap[i,j,k]
            end
            TKE[i,j,k] = 0.5*(u[i,j,k]^2 + v[i,j,k]^2 + w[i,j,k]^2)
        end
    end
end
mu_eff_2 = sum(mu_eff[:,:,:])/(n^3)

# Compute spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(mu_eff))/(N^3)
Ek = scalar_hat.^2
Ek = fftshift(Ek)

box_sidex = N
box_sidey = N
box_sidez = N

box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

centerx = Int(box_sidex/2)
centery = Int(box_sidey/2)
centerz = Int(box_sidez/2)

Ek_avsphr = zeros(box_radius-1)
wavenumber = zeros(box_radius-1)

for k = 1:N
	for j = 1:N
		for i = 1:N
            global m
			m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
            if (m <= box_radius-1)
                Ek_avsphr[m] = Ek_avsphr[m] + Ek[i,j,k]
                wavenumber[m] = (2*pi*m/L)
            end
        end
    end
end

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :blue4,linestyle =:dash,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,legend=false,xticks=1:0.5:3.5,xlims=[1,3.5])


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var3 = []
var4 = []
var5 = []
var6 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_H_1/nodal_values/nodal_values-25.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,5])
    append!(var5,VAR[:,7])
    append!(var6,VAR[:,8])
end

n = 192
u = zeros(n,n,n)
v = zeros(n,n,n)
w = zeros(n,n,n)
scalar = zeros(n,n,n)
forcing = zeros(n,n,n)
latent = zeros(n,n,n)

(u,count) = nodal_value(var1,u,Nx,Ny,Nz)
(v,count) = nodal_value(var2,v,Nx,Ny,Nz)
(w,count) = nodal_value(var3,w,Nx,Ny,Nz)
(scalar,count) = nodal_value(var4,scalar,Nx,Ny,Nz)
(forcing,count) = nodal_value(var5,forcing,Nx,Ny,Nz)
(latent,count) = nodal_value(var6,latent,Nx,Ny,Nz)


# Calculate the advection term
L = 0.512
h = L/(n-1)

advec = zeros(n,n,n)
for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            global advec
            advec[i,j,k] = u[i,j,k]*((scalar[i+1,j,k] - scalar[i,j,k])/h) + v[i,j,k]*((scalar[i,j+1,k] - scalar[i,j,k])/h) + w[i,j,k]*((scalar[i,j,k+1] - scalar[i,j,k])/h)
        end
    end
end
advec[n,n,n] = u[n,n,n]*((scalar[1,n,n] - scalar[n,n,n])/h) + v[n,n,n]*((scalar[n,1,n] - scalar[n,n,n])/h) + w[n,n,n]*((scalar[n,n,1] - scalar[n,n,n])/h)

# Calculate the Laplacian of the scalar
lap = zeros(n,n,n)
for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global lap
            lap[i-1,j-1,k-1] = (scalar[i+1,j,k] - 2*scalar[i,j,k] + scalar[i-1,j,k])/(h^2) + (scalar[i,j+1,k] - 2*scalar[i,j,k] + scalar[i,j-1,k])/(h^2) + (scalar[i,j,k+1] - 2*scalar[i,j,k] + scalar[i,j,k-1])/(h^2)
        end
    end
end

lap[1,1,1] = (scalar[2,1,1] - 2*scalar[1,1,1] + scalar[n,1,1])/(h^2) + (scalar[1,2,1] - 2*scalar[1,1,1] + scalar[1,n,1])/(h^2) + (scalar[1,1,2] - 2*scalar[1,1,1] + scalar[1,1,n])/(h^2)

lap[n,n,n] = (scalar[1,n,n] - 2*scalar[n,n,n] + scalar[n-1,n,n])/(h^2) + (scalar[n,1,n] - 2*scalar[n,n,n] + scalar[n,n-1,n])/(h^2) + (scalar[n,n,1] - 2*scalar[n,n,n] + scalar[n,n,n-1])/(h^2)


# Calculate artificial forcing viscosity, TKE, and scalar variance
mu_eff = zeros(n,n,n)
TKE = zeros(n,n,n)
mu_v = 2.2e-5
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff, TKE, count
            if (lap[i,j,k] != 0)
                mu_eff[i,j,k] = (advec[i,j,k] - latent[i,j,k] - forcing[i,j,k])/lap[i,j,k]
            end
            TKE[i,j,k] = 0.5*(u[i,j,k]^2 + v[i,j,k]^2 + w[i,j,k]^2)
        end
    end
end
mu_eff_3 = sum(mu_eff[:,:,:])/(n^3)
#=
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

mu_eff_diag = zeros(n)
TKE_diag = zeros(n)
forcing_diag = zeros(n)
latent_diag = zeros(n)
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_forcing_diag, count
            if (i == j && j == k)
                count = count + 1
                mu_eff_diag[count] = mu_eff[i,j,k] 
                TKE_diag[count] = TKE[i,j,k] 
                forcing_diag[count] = forcing[i,j,k]
                latent_diag[count] = latent[i,j,k]  
            end
        end
    end
end

X = collect(0:c:D)

plot(X[1:end],mu_eff_diag[1:end],color = :green4,linestyle =:dashdot,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,legend=false)
=#

# Compute spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(mu_eff))/(N^3)
Ek = scalar_hat.^2
Ek = fftshift(Ek)

box_sidex = N
box_sidey = N
box_sidez = N

box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

centerx = Int(box_sidex/2)
centery = Int(box_sidey/2)
centerz = Int(box_sidez/2)

Ek_avsphr = zeros(box_radius-1)
wavenumber = zeros(box_radius-1)

for k = 1:N
	for j = 1:N
		for i = 1:N
            global m
			m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
            if (m <= box_radius-1)
                Ek_avsphr[m] = Ek_avsphr[m] + Ek[i,j,k]
                # wavenumber[m] = (2*pi*m/L)*(eta/sqrt(0.68))
                # wavenumber[m] = (2*pi*m/L)*eta
                wavenumber[m] = (2*pi*m/L)
            end
        end
    end
end

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :green4,linestyle =:dashdot,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var3 = []
var4 = []
var5 = []
var6 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_1/nodal_values/nodal_values-25.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,5])
    append!(var5,VAR[:,7])
    append!(var6,VAR[:,8])
end

n = 192
u = zeros(n,n,n)
v = zeros(n,n,n)
w = zeros(n,n,n)
scalar = zeros(n,n,n)
forcing = zeros(n,n,n)
latent = zeros(n,n,n)

(u,count) = nodal_value(var1,u,Nx,Ny,Nz)
(v,count) = nodal_value(var2,v,Nx,Ny,Nz)
(w,count) = nodal_value(var3,w,Nx,Ny,Nz)
(scalar,count) = nodal_value(var4,scalar,Nx,Ny,Nz)
(forcing,count) = nodal_value(var5,forcing,Nx,Ny,Nz)
(latent,count) = nodal_value(var6,latent,Nx,Ny,Nz)


# Calculate the advection term
L = 0.512
h = L/(n-1)

advec = zeros(n,n,n)
for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            global advec
            advec[i,j,k] = u[i,j,k]*((scalar[i+1,j,k] - scalar[i,j,k])/h) + v[i,j,k]*((scalar[i,j+1,k] - scalar[i,j,k])/h) + w[i,j,k]*((scalar[i,j,k+1] - scalar[i,j,k])/h)
        end
    end
end
advec[n,n,n] = u[n,n,n]*((scalar[1,n,n] - scalar[n,n,n])/h) + v[n,n,n]*((scalar[n,1,n] - scalar[n,n,n])/h) + w[n,n,n]*((scalar[n,n,1] - scalar[n,n,n])/h)


# Calculate the Laplacian of the scalar
lap = zeros(n,n,n)
for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global lap
            lap[i-1,j-1,k-1] = (scalar[i+1,j,k] - 2*scalar[i,j,k] + scalar[i-1,j,k])/(h^2) + (scalar[i,j+1,k] - 2*scalar[i,j,k] + scalar[i,j-1,k])/(h^2) + (scalar[i,j,k+1] - 2*scalar[i,j,k] + scalar[i,j,k-1])/(h^2)
        end
    end
end

lap[1,1,1] = (scalar[2,1,1] - 2*scalar[1,1,1] + scalar[n,1,1])/(h^2) + (scalar[1,2,1] - 2*scalar[1,1,1] + scalar[1,n,1])/(h^2) + (scalar[1,1,2] - 2*scalar[1,1,1] + scalar[1,1,n])/(h^2)

lap[n,n,n] = (scalar[1,n,n] - 2*scalar[n,n,n] + scalar[n-1,n,n])/(h^2) + (scalar[n,1,n] - 2*scalar[n,n,n] + scalar[n,n-1,n])/(h^2) + (scalar[n,n,1] - 2*scalar[n,n,n] + scalar[n,n,n-1])/(h^2)


# Calculate artificial forcing viscosity, TKE, and scalar variance
mu_eff = zeros(n,n,n)
mu_forcing = zeros(n,n,n)
TKE = zeros(n,n,n)
mu_v = 2.2e-5
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff, mu_forcing, TKE, count
            if (lap[i,j,k] != 0)
                mu_eff[i,j,k] = (advec[i,j,k] - latent[i,j,k] - forcing[i,j,k])/lap[i,j,k]
            end
            TKE[i,j,k] = 0.5*(u[i,j,k]^2 + v[i,j,k]^2 + w[i,j,k]^2)
        end
    end
end
mu_eff_4 = sum(mu_eff[:,:,:])/(n^3)

#=
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

mu_eff_diag = zeros(n)
TKE_diag = zeros(n)
forcing_diag = zeros(n)
latent_diag = zeros(n)
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_forcing_diag, count
            if (i == j && j == k)
                count = count + 1
                mu_eff_diag[count] = mu_eff[i,j,k] 
                TKE_diag[count] = TKE[i,j,k] 
                forcing_diag[count] = forcing[i,j,k]
                latent_diag[count] = latent[i,j,k]  
            end
        end
    end
end

X = collect(0:c:D)

plot(X[1:end],mu_eff_diag[1:end],color = :orange3,linestyle =:solid,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)
=#

# Compute spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(mu_eff))/(N^3)
Ek = scalar_hat.^2
Ek = fftshift(Ek)

box_sidex = N
box_sidey = N
box_sidez = N

box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

centerx = Int(box_sidex/2)
centery = Int(box_sidey/2)
centerz = Int(box_sidez/2)

Ek_avsphr = zeros(box_radius-1)
wavenumber = zeros(box_radius-1)

for k = 1:N
	for j = 1:N
		for i = 1:N
            global m
			m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
            if (m <= box_radius-1)
                Ek_avsphr[m] = Ek_avsphr[m] + Ek[i,j,k]
                # wavenumber[m] = (2*pi*m/L)*(eta/sqrt(0.68))
                # wavenumber[m] = (2*pi*m/L)*eta
                wavenumber[m] = (2*pi*m/L)
            end
        end
    end
end

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :orange3,linestyle =:solid,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)


# Read data
N = 128
Nx = 8
Ny = 4
Nz = 4

var1 = []
var2 = []
var3 = []
var4 = []
var5 = []
var6 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_L_1/nodal_values/nodal_values-25.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
    append!(var4,VAR[:,5])
    append!(var5,VAR[:,7])
    append!(var6,VAR[:,8])
end

n = 192
u = zeros(n,n,n)
v = zeros(n,n,n)
w = zeros(n,n,n)
scalar = zeros(n,n,n)
forcing = zeros(n,n,n)
latent = zeros(n,n,n)

(u,count) = nodal_value(var1,u,Nx,Ny,Nz)
(v,count) = nodal_value(var2,v,Nx,Ny,Nz)
(w,count) = nodal_value(var3,w,Nx,Ny,Nz)
(scalar,count) = nodal_value(var4,scalar,Nx,Ny,Nz)
(forcing,count) = nodal_value(var5,forcing,Nx,Ny,Nz)
(latent,count) = nodal_value(var6,latent,Nx,Ny,Nz)


# Calculate the advection term
L = 0.512
h = L/(n-1)

advec = zeros(n,n,n)
for k = 1:n-1
    for j = 1:n-1
        for i = 1:n-1
            global advec
            advec[i,j,k] = u[i,j,k]*((scalar[i+1,j,k] - scalar[i,j,k])/h) + v[i,j,k]*((scalar[i,j+1,k] - scalar[i,j,k])/h) + w[i,j,k]*((scalar[i,j,k+1] - scalar[i,j,k])/h)
        end
    end
end
advec[n,n,n] = u[n,n,n]*((scalar[1,n,n] - scalar[n,n,n])/h) + v[n,n,n]*((scalar[n,1,n] - scalar[n,n,n])/h) + w[n,n,n]*((scalar[n,n,1] - scalar[n,n,n])/h)


# Calculate the Laplacian of the scalar
lap = zeros(n,n,n)
for k = 2:n-1
    for j = 2:n-1
        for i = 2:n-1
            global lap
            lap[i-1,j-1,k-1] = (scalar[i+1,j,k] - 2*scalar[i,j,k] + scalar[i-1,j,k])/(h^2) + (scalar[i,j+1,k] - 2*scalar[i,j,k] + scalar[i,j-1,k])/(h^2) + (scalar[i,j,k+1] - 2*scalar[i,j,k] + scalar[i,j,k-1])/(h^2)
        end
    end
end

lap[1,1,1] = (scalar[2,1,1] - 2*scalar[1,1,1] + scalar[n,1,1])/(h^2) + (scalar[1,2,1] - 2*scalar[1,1,1] + scalar[1,n,1])/(h^2) + (scalar[1,1,2] - 2*scalar[1,1,1] + scalar[1,1,n])/(h^2)

lap[n,n,n] = (scalar[1,n,n] - 2*scalar[n,n,n] + scalar[n-1,n,n])/(h^2) + (scalar[n,1,n] - 2*scalar[n,n,n] + scalar[n,n-1,n])/(h^2) + (scalar[n,n,1] - 2*scalar[n,n,n] + scalar[n,n,n-1])/(h^2)


# Calculate artificial forcing viscosity, TKE, and scalar variance
mu_eff = zeros(n,n,n)
mu_forcing = zeros(n,n,n)
TKE = zeros(n,n,n)
mu_v = 2.2e-5
count = 0
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff, mu_forcing, TKE, count
            if (lap[i,j,k] != 0)
                mu_eff[i,j,k] = (advec[i,j,k] - latent[i,j,k] - forcing[i,j,k])/lap[i,j,k]
            end
            TKE[i,j,k] = 0.5*(u[i,j,k]^2 + v[i,j,k]^2 + w[i,j,k]^2)
        end
    end
end
mu_eff_5 = sum(mu_eff[:,:,:])/(n^3)
#=
D = sqrt(L^2+L^2+L^2)
c = D/(n-1)

mu_eff_diag = zeros(n)
TKE_diag = zeros(n)
forcing_diag = zeros(n)
latent_diag = zeros(n)
for k = 1:n
    for j = 1:n
        for i = 1:n
            global mu_eff_diag, count
            if (i == j && j == k)
                count = count + 1
                mu_eff_diag[count] = mu_eff[i,j,k] 
                TKE_diag[count] = TKE[i,j,k]
                forcing_diag[count] = forcing[i,j,k] 
                latent_diag[count] = latent[i,j,k]  
            end
        end
    end
end

X = collect(0:c:D)

plot!(X[1:end],mu_eff_diag[1:end],color = :black,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)
scatter!(X[1:6:end],mu_eff_diag[1:6:end],color = :black,markershape = :square,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)
=#

# Compute spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(mu_eff))/(N^3)
Ek = scalar_hat.^2
Ek = fftshift(Ek)

box_sidex = N
box_sidey = N
box_sidez = N

box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

centerx = Int(box_sidex/2)
centery = Int(box_sidey/2)
centerz = Int(box_sidez/2)

Ek_avsphr = zeros(box_radius-1)
wavenumber = zeros(box_radius-1)

for k = 1:N
	for j = 1:N
		for i = 1:N
            global m
			m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
            if (m <= box_radius-1)
                Ek_avsphr[m] = Ek_avsphr[m] + Ek[i,j,k]
                # wavenumber[m] = (2*pi*m/L)*(eta/sqrt(0.68))
                # wavenumber[m] = (2*pi*m/L)*eta
                wavenumber[m] = (2*pi*m/L)
            end
        end
    end
end

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :black,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5)

scatter!(log.(10,wavenumber[1:6:end]),log.(10,Ek_avsphr[1:6:end]),color = :black,markershape = :square,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)

