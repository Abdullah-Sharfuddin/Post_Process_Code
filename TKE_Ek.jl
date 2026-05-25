using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using FFTW


function nodal_value(var,vel,Nx,Ny,Nz)

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
                            vel[i,j,k] = var[count]
                        end
                    end
                end
            end
        end
    end

    return vel
end


function rms_spectra(vel_x,vel_y,vel_z)

    # Compute fluctuation
    mean_vel_x = sum(vel_x[:,:,:]) / (n*n*n)
    mean_vel_y = sum(vel_y[:,:,:]) / (n*n*n)
    mean_vel_z = sum(vel_z[:,:,:]) / (n*n*n)
    for k = 1:n
        for j = 1:n
            for i = 1:n
                vel_x[i,j,k] = vel_x[i,j,k] - mean_vel_x
                vel_y[i,j,k] = vel_y[i,j,k] - mean_vel_y
                vel_z[i,j,k] = vel_z[i,j,k] - mean_vel_z
            end
        end
    end

    # Compute rms fluctuation of combined velocity field
    ms_vel_x = sum(vel_x[:,:,:].^2) / (n*n*n)
    ms_vel_y = sum(vel_y[:,:,:].^2) / (n*n*n)
    ms_vel_z = sum(vel_z[:,:,:].^2) / (n*n*n)

    Urms = sqrt(ms_vel_x + ms_vel_y + ms_vel_z)

    # Compute TKE in physical space
    E = zeros(n,n,n)
    for k = 1:n
        for j = 1:n
            for i = 1:n
                E[i,j,k] = 0.5*(vel_x[i,j,k]^2 + vel_y[i,j,k]^2 + vel_z[i,j,k]^2)
            end
        end
    end

    # Compute TKE in Fourier space
    N = n
    L = 0.512

    vel_x_hat = abs.(fft(vel_x))/(N^3)
    vel_y_hat = abs.(fft(vel_y))/(N^3)
    vel_z_hat = abs.(fft(vel_z))/(N^3)

    Ek_U = vel_x_hat.^2
    Ek_V = vel_y_hat.^2
    Ek_W = vel_z_hat.^2

    Ek_U = fftshift(Ek_U)
    Ek_V = fftshift(Ek_V)
    Ek_W = fftshift(Ek_W)

    box_sidex = N
    box_sidey = N
    box_sidez = N

    box_radius = Int(ceil((sqrt((box_sidex)^2+(box_sidey)^2+(box_sidez)^2))/2)+1)

    centerx = Int(box_sidex/2)
    centery = Int(box_sidey/2)
    centerz = Int(box_sidez/2)

    Ek_U_avsphr = zeros(box_radius-1)
    Ek_V_avsphr = zeros(box_radius-1)
    Ek_W_avsphr = zeros(box_radius-1)

    wavenumber = zeros(box_radius-1)

    for k = 1:N
        for j = 1:N
            for i = 1:N
                global m
                m =  Int(round(sqrt((i-centerx)^2+(j-centery)^2+(k-centerz)^2))+1)
                if (m <= box_radius-1)
                    Ek_U_avsphr[m] = Ek_U_avsphr[m] + Ek_U[i,j,k]
                    Ek_V_avsphr[m] = Ek_V_avsphr[m] + Ek_V[i,j,k]
                    Ek_W_avsphr[m] = Ek_W_avsphr[m] + Ek_W[i,j,k]
                    wavenumber[m] = 2*pi*m/L
                end
            end
        end
    end

    Ek_avsphr = 0.5*(Ek_U_avsphr + Ek_V_avsphr + Ek_W_avsphr)

    return Urms, wavenumber, Ek_avsphr
end

#=
# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseMC1_256/nodal_values/nodal_values-5.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms1,wavenumber1,Ek_avsphr1) = rms_spectra(vel_x,vel_y,vel_z)

eta1 = 0.00147

nu = 1.5e-5
wavenumber1 = wavenumber1 .* eta1
# Ek_avsphr1 = Ek_avsphr1 ./ ((nu^(5/4)).*(eps1^(1/4)))

plot(log.(10,wavenumber1),log.(10,Ek_avsphr1),color = :red4,linestyle =:dot,linewidth=2.5,
labels="L",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))

# X = log.(10,wavenumber1)
# Y = log.(10,Ek_avsphr1)

# m = (Y[80] - Y[20]) / (X[80] - X[20])
# b = Y[20] - m*X[20]
# plot!(x -> m*x + b, X[20], X[80],color = :purple,linestyle =:dot,linewidth=2.5,
# labels="-17/3",grid=false,thickness_scaling=1.5,legend=false)



# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/nodal_values/nodal_values-5.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms2,wavenumber2,Ek_avsphr2) = rms_spectra(vel_x,vel_y,vel_z)

eta2 = 0.00103

# eta2 = 0.00169
# eps2 = 0.00041

nu = 1.5e-5
wavenumber2 = wavenumber2 .* eta2
# Ek_avsphr2 = Ek_avsphr2 ./ ((nu^(5/4)).*(eps2^(1/4)))

plot!(log.(10,wavenumber2),log.(10,Ek_avsphr2),color = :blue4,linestyle =:dash,linewidth=2,
labels="M",grid=false,thickness_scaling=1.5)
=#

#=
# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseHC1_256/nodal_values/nodal_values-0.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms3,wavenumber3,Ek_avsphr3) = rms_spectra(vel_x,vel_y,vel_z)

eta3 = 0.000745
eps3 = 0.010927

# eta3 = 0.00105
# eps3 = 0.00273

nu = 1.5e-5
wavenumber3 = wavenumber3 .* eta3
# Ek_avsphr3 = Ek_avsphr3 ./ ((nu^(5/4)).*(eps3^(1/4)))

plot!(log.(10,wavenumber3),log.(10,Ek_avsphr3),color = :green4,linestyle =:solid,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,legend=false)


# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseVC1_256/nodal_values/nodal_values-0.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms4,wavenumber4,Ek_avsphr4) = rms_spectra(vel_x,vel_y,vel_z)

eta4 = 0.000515
eps4 = 0.047908

# eta4 = 0.00073
# eps4 = 0.00118

nu = 1.5e-5
wavenumber4 = wavenumber4 .* eta4
# Ek_avsphr4 = Ek_avsphr4 ./ ((nu^(5/4)).*(eps4^(1/4)))

plot!(log.(10,wavenumber4),log.(10,Ek_avsphr4),color = :black,linestyle =:dashdot,linewidth=2,ylims=[-20,0],
labels="V",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
=#


#=
val1 = readdlm("Post_Process_Scalar/Data/M_TKE.txt",skipstart=1)
plot!(val1[:,1],val1[:,2],color = :blue4,linestyle =:dash,linewidth=2,
labels="M",grid=false,thickness_scaling=1.5)
=#


#=
val3 = readdlm("Post_Process_Scalar/Data/64_TKE.txt",skipstart=1)
plot(val3[:,1],val3[:,2] .- 6.5,color = :nuuk,linestyle =:dot,linewidth=2.5,
labels="64",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))

val2 = readdlm("Post_Process_Scalar/Data/128_TKE.txt",skipstart=1)
plot!(val2[:,1],val2[:,2] .- 6.5,color = :purple,linestyle =:dashdot,linewidth=2,
labels="128",grid=false,thickness_scaling=1.5)


# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Scalar/F3_CaseMC1_256/nodal_values/nodal_values-5.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms,wavenumber,Ek_avsphr) = rms_spectra(vel_x,vel_y,vel_z)

eta = 0.00136
eps = 0.000986

nu = 1.5e-5
wavenumber = wavenumber .* eta
Ek_avsphr = Ek_avsphr ./ ((nu^(5/4)).*(eps^(1/4)))

plot!(log.(10,wavenumber),log.(10,Ek_avsphr) .- 6.5,color = :summer,linestyle =:dash,linewidth=2,
labels="256",grid=false,thickness_scaling=1.5,legend=false)

val1 = readdlm("Post_Process_Scalar/Data/512_TKE.txt",skipstart=1)
plot!(val1[:,1],val1[:,2] .- 6.5,color = :black,linestyle =:solid,linewidth=2,
labels="512",grid=false,thickness_scaling=1.5,yticks=-15:2.5:-2.5,ylims=[-15,-2.5])
=#




# Read data
N = 1024
Nx = 16
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Upscale/Case_H_P_256/nodal_values/nodal_values-10.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms4,wavenumber4,Ek_avsphr4) = rms_spectra(vel_x,vel_y,vel_z)

# eta4 = 0.00073
# eps4 = 0.00118

# nu = 1.5e-5
# wavenumber4 = wavenumber4 .* eta4

plot(log.(10,wavenumber4),log.(10,Ek_avsphr4),color = :summer,linestyle =:solid,linewidth=2,
labels="H",grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legend=false,xticks=1:0.5:3.5,xlims=[1,3.6])


#=
# Read data
N = 512
Nx = 8
Ny = 8
Nz = 8

var1 = []
var2 = []
var3 = []
I = 0
for i = 1:N
    global I, var
    I = i - 1
    VAR = readdlm("PR_DNS_Cloud/F3_CaseD_128/nodal_values/nodal_values-6.00-$I")
    append!(var1,VAR[:,1])
    append!(var2,VAR[:,2])
    append!(var3,VAR[:,3])
end

n = 256
vel_x = zeros(n,n,n)
vel_y = zeros(n,n,n)
vel_z = zeros(n,n,n)

vel_x = nodal_value(var1,vel_x,Nx,Ny,Nz)
vel_y = nodal_value(var2,vel_y,Nx,Ny,Nz)
vel_z = nodal_value(var3,vel_z,Nx,Ny,Nz)

(Urms4,wavenumber4,Ek_avsphr4) = rms_spectra(vel_x,vel_y,vel_z)

eta4 = 0.0020

# eta4 = 0.00073
# eps4 = 0.00118

nu = 1.5e-5
wavenumber4 = wavenumber4 .* eta4

plot!(log.(10,wavenumber4),log.(10,Ek_avsphr4),color = :bwr,linestyle =:dashdot,linewidth=2,
labels=false,grid=false,thickness_scaling=1.5,tickfont=font(9,"Helvetica Bold"),legendfont=font(5,"Helvetica Bold"))
=#


