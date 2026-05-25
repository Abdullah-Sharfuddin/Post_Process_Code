using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW

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
N = 128
Nx = 8
Ny = 4
Nz = 4

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.0019

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

plot(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :red4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))
# xticks=-2:0.4:0.6,xlims=[-2,0.6]
# Match index is 64 for H cases and 30 for L cases


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
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :blue4,linestyle =:dash,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,legend=false,xticks=1:0.5:3.5,xlims=[1,3.5])



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
    VAR = readdlm("PR_DNS_Scalar/Case_L_H_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/Case_H_L_2/nodal_values/nodal_values-25.00-$I")
    append!(var,VAR[:,6])
end

n = 192
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

# plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :black,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,xticks=1:0.5:3.5,xlims=[1,3.5])
scatter!(log.(10,wavenumber[1:6:end]),log.(10,Ek_avsphr[1:6:end]),color = :black,markershape = :square,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)
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

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

plot(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :red4,linestyle =:dot,linewidth=2.5,
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)



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
    VAR = readdlm("PR_DNS_Mixed/Case_H_S_H/nodal_values/nodal_values-30.00-$I")
    append!(var,VAR[:,2])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

# plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :orange3,linestyle =:solid,linewidth=2,
# labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :bwr,linestyle =:dash,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,xticks=1:0.5:3.5,xlims=[1,3.5])
scatter!(log.(10,wavenumber[1:6:end]),log.(10,Ek_avsphr[1:6:end]),color = :bwr,markershape = :circle,markersize = :2,
labels="Critical",grid=false,thickness_scaling=1.5)

# A = readdlm("Post_Process_Scalar/Data/LA_L_1.txt",skipstart=1)
# plot!(A[:,1],A[:,2],color = :yellow4,linestyle =:solid,linewidth=2,yticks=-20:5:-5,ylims=[-20,-5],
# labels="V",grid=false,thickness_scaling=1.5,legendfont=font(5,"Helvetica Bold"))
=#

#=
# Thermal Boundary
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
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_128/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,9])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.0019

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

plot(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :red4,linestyle =:dot,linewidth=2.5,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"))


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
    VAR = readdlm("PR_DNS_Upscale/Case_H_D_128/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,9])
end

n = 128
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
# eta = 0.00089

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

plot!(log.(10,wavenumber[1:end]),log.(10,Ek_avsphr[1:end]),color = :blue4,linestyle =:dash,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,legend=false)
=#



## Upscale # Spectral Break
N = 1024
Nx = 16
Ny = 8
Nz = 8

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,9])
end

n = 256
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# Compute liquid water mixing ratio
rho_l = 1000
h = 0.512/256
rho_a = 1.0
C = (4*pi*rho_l)/(3*rho_a*(h^3))
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = C*(scalar[i,j,k]^3)
        end
    end
end


# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
eta = 0.00088196095823
x = log.(10,wavenumber[1:end] .*eta)
y = log.(10,Ek_avsphr[1:end])

x_ticks_pre = -2:0.5:0.5
x_labels = ["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}", "10^{0.5}"]

# y_ticks_pre = -15:3:-3
# y_labels = ["10^{-15}", "10^{-12}", "10^{-9}", "10^{-6}", "10^{-3}"]

y_ticks_pre = -13:1:-9
y_labels = ["10^{-13}", "10^{-12}", "10^{-11}", "10^{-10}", "10^{-9}"]

plot(x,y,color = :summer,linestyle =:solid,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)



# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var,VAR[:,9])
end

n = 256
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# Compute liquid water mixing ratio
rho_l = 1000
h = 0.512/256
rho_a = 1.0
C = (4*pi*rho_l)/(3*rho_a*(h^3))
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = C*(scalar[i,j,k]^3)
        end
    end
end


# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
eta = 0.00153371925558
x = log.(10,wavenumber[1:end] .*eta)
y = log.(10,Ek_avsphr[1:end])

x_ticks_pre = -2:0.5:0.5
x_labels = ["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}", "10^{0.5}"]

# y_ticks_pre = -15:3:-3
# y_labels = ["10^{-15}", "10^{-12}", "10^{-9}", "10^{-6}", "10^{-3}"]

y_ticks_pre = -13:1:-9
y_labels = ["10^{-13}", "10^{-12}", "10^{-11}", "10^{-10}", "10^{-9}"]

plot!(x,y,color = :bwr,linestyle =:dash,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)


# x1 = zeros(222); y1 = zeros(222)
# y_ticks_pre = -18:1:-13
# y_labels = ["10^{-18}", "10^{-17}", "10^{-16}", "10^{-15}", "10^{-14}", "10^{-13}"]

# plot!(twinx(),x1,y1,color = :bwr,linestyle =:dash,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
# labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),ylims = (-18, -13))


#=
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
    VAR = readdlm("PR_DNS_Upscale/Case_H_P0_256/nodal_values/nodal_values-8.00-$I")
    append!(var,VAR[:,8])
end

n = 256
scalar = zeros(n,n,n)

(scalar,count) = nodal_value(var,scalar,Nx,Ny,Nz)


# # Compute liquid water mixing ratio
# rho_l = 1000
# h = 0.512/256
# rho_a = 1.0
# C = (4*pi*rho_l)/(3*rho_a*(h^3))
# for k = 1:n
#     for j = 1:n
#         for i = 1:n
#             scalar[i,j,k] = C*(scalar[i,j,k]^3)
#         end
#     end
# end


# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
eta = 0.00088567694246
x = log.(10,wavenumber[1:end] .*eta)
y = log.(10,Ek_avsphr[1:end])

x_ticks_pre = -2:0.5:0.5
x_labels = ["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}", "10^{0.5}"]

# y_ticks_pre = -15:3:-3
# y_labels = ["10^{-15}", "10^{-12}", "10^{-9}", "10^{-6}", "10^{-3}"]

y_ticks_pre = -12:2:-4
y_labels = ["10^{-12}", "10^{-10}", "10^{-8}", "10^{-6}", "10^{-4}"]

plot!(x,y,color = :purple,linestyle =:dot,linewidth=2.5,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)
=#

#=
# Spectra of particle radii
# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
I = 0
for i = 1:N
    global I, var1
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/particle/particle_radius-10.00-$I")
    append!(var1,VAR[:,2])
end

n = 250
scalar = zeros(n,n,n)

# Assign nodal values
count = 1
for k = 1:n
    for j = 1:n
        for i = 1:n
            global count
            scalar[i,j,k] = var1[count]
            count = count + 1
        end
    end
end

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
eta = 0.00088196095823
x = log.(10,wavenumber[1:end] .*eta)
y = log.(10,Ek_avsphr[1:end])

x_ticks_pre = -2:0.5:0.5
x_labels = ["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}", "10^{0.5}"]

y_ticks_pre = -21:2:-13
y_labels = ["10^{-21}", "10^{-19}", "10^{-17}", "10^{-15}", "10^{-13}"]

plot(x,y,color = :summer,linestyle =:solid,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)


# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
I = 0
for i = 1:N
    global I, var1
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_L_P_256/particle/particle_radius-10.00-$I")
    append!(var1,VAR[:,2])
end

n = 250
scalar = zeros(n,n,n)

# Assign nodal values
count = 1
for k = 1:n
    for j = 1:n
        for i = 1:n
            global count
            scalar[i,j,k] = var1[count]
            count = count + 1
        end
    end
end

# Compute mean scalar field
mean_scalar = sum(scalar[:,:,:]) / (n*n*n)

# Compute fluctuating scalar field
for k = 1:n
    for j = 1:n
        for i = 1:n
            scalar[i,j,k] = scalar[i,j,k] - mean_scalar
        end
    end
end

# Compute scalar spectrum
N = n
L = 0.512

scalar_hat = abs.(fft(scalar))/(N^3)
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
eta = 0.00153371925558
x = log.(10,wavenumber[1:end] .*eta)
y = log.(10,Ek_avsphr[1:end])

x_ticks_pre = -2:0.5:0.5
x_labels = ["10^{-2}", "10^{-1.5}", "10^{-1}", "10^{-0.5}", "10^{0}", "10^{0.5}"]

y_ticks_pre = -21:2:-13
y_labels = ["10^{-21}", "10^{-19}", "10^{-17}", "10^{-15}", "10^{-13}"]

plot!(x,y,color = :bwr,linestyle =:dash,linewidth=2,xticks = (x_ticks_pre, x_labels), yticks = (y_ticks_pre, y_labels),
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false)
=#