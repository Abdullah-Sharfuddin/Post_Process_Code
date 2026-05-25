using DelimitedFiles
using Distributions
using Plots
using LaTeXStrings
using LinearAlgebra
using SparseArrays

# Variational approach
# 3D VAR implementation

# # Arguments
# x_background : Background state vector (n × 1)
# y_obs : Observation vector (p × 1)
# H: Observation operator matrix (p × n) or function H(x) that maps state to observation space
# R: Observation error covariance matrix (p × p)
# B: Background error covariance matrix (n × n)
# method : Solution method (:direct, :cg for conjugate gradient)

# # Returns
# x_analysis : Analysis state vector (n × 1)


function var3d(x_background, y_obs, H, R, B; method=:direct)
    # Convert H to matrix if it's a function
    if isa(H, Function)
        H_matrix = calculate_jacobian(H, x_background)
    else
        H_matrix = H
    end
    
    # Difference vector
    d = y_obs - H_matrix * x_background
    
    # Solve for x_analysis by minimizing the cost function
    if method == :direct
        # Direct solution using matrix inversion
        BHT = B * H_matrix'
        HBHT = H_matrix * BHT
        # Compute Kalman gain matrix
        K = BHT / (HBHT + R)
        # Compute analysis vector
        x_analysis = x_background + K * d

    elseif method == :cg
        # Iterative solution using conjugate gradient
        # Solve (H B H' + R) w = d for w
        A = H_matrix * B * H_matrix' + R
        w = conjugate_gradient(A, d)
        # Compute analysis vector
        x_analysis = x_background + B * H_matrix' * w
    else
        error("Method not specified: $method. Use :direct or :cg.")
    end
    
    return x_analysis
end


# Numerically calculate Jacobian matrix of observation operator H at point x.
function calculate_jacobian(H, x)
    n = length(x)
    p = length(H(x))
    J = zeros(p, n)
    ϵ = 1e-6
    
    for i in 1:n
        dx = zeros(n)
        dx[i] = ϵ
        J[:, i] = (H(x + dx) - H(x - dx)) / (2ϵ)
    end
    
    return J
end

# Conjugate gradient method for solving Ax = b.
function conjugate_gradient(A, b)
    n = size(b,1)
    x0 = zeros(n)
    r = b - A * x0
    p = r
    rsold = r' * r
    x = x0

    max_iter=1000
    tol=1e-6
    
    while sqrt(rsold) > tol
        Ap = A * p
        alpha = rsold / (p' * Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = r' * r
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    end
    
    return x
end


# Create the B matrix with exponential decay based on distance.
# Gaussian correlation model
# The B matrix is usually non-diagonal due to spatial and multivariate correlations
function B_matrix(n, length_scale, variance)
    # For simplicity, assume a vector
    grid_points = collect(1:n)
    B = zeros(n, n)
    
    for i in 1:n
        for j in 1:n
            distance = abs(grid_points[i] - grid_points[j])
            B[i,j] = variance * exp(-distance / length_scale)
        end
    end
    
    return B
end


# Number of elements
n = 20000  # For background state
p = 2000   # For external observation

# Create background state (Should come from numerical model)
x_true = randn(n)
x_background = x_true + randn(n) * 0.5  # Add some error

# Observation operator (simple linear operator for this example)
H = sprandn(p, n, 0.2)  # Sparse random matrix

# Generate observations (Should come from external observation: radar, satellite etc.)
y_true = H * x_true
y_obs = y_true + randn(p) * 0.2  # Add some error

# Create covariance matrices
# The R matrix represents the uncertainties associated with observations
# Assumed to be diagonal because observation errors are uncorrelated at different observation points
# Observation error covariance (0.2^2)
R = Diagonal(ones(p) * 0.04)

# The B matrix captures how errors in one model variable at a location influence errors at other locations
# Background error covariance (0.5^2)
B = B_matrix(n, 10.0, 0.25)

# Perform 3D-VAR
x_analysis_direct = var3d(x_background, y_obs, H, R, B, method=:direct)
x_analysis_cg = var3d(x_background, y_obs, H, R, B, method=:cg)

# Calculate errors
background_error = norm(x_background - x_true)
analysis_error_direct = norm(x_analysis_direct - x_true)
analysis_error_cg = norm(x_analysis_cg - x_true)

## Replace x_background with numerical data for temperature/humidity/particle radius
## Replace y_obs with appropriate radar/satellite data when available

