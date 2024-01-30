using LinearAlgebra


#---------------------------------------#
#               Outline                 #
#---------------------------------------#

# - Parameters
# - Differential Equation
# - Quadric RBFs & Derivatives
# - Iterator Points
# - Linear System Generator
# - Approximate Solution 
# - Example Use


#---------------------------------------#
#              Parameters               #
#---------------------------------------#

Δ = 3      # Number of partitions along one-dimension, the total number of 
            # points along a dimension is then Δ + 1 and includes the boundary.
            # (Δ+1)^3 is the total number of grid points in three dimensions.
            # Δ = 5000, should be 932 Gigabytes

σ = 1.0     # Decay factor for the quadric radial basis functions 


#---------------------------------------#
#         Differential Equation         #
#---------------------------------------#

# Oracle Solution 
u_sol(ξ) = sin(π*ξ[1]) * sin(π*ξ[2]/2) * sin(3*π*ξ[3]/2)

# We study Poisson's equation, Δu(ξ) = h(ξ) on the unit cube
h(ξ) = (-7*π^2/2) * sin(π*ξ[1]) * sin(π*ξ[2]/2) * sin(3*π*ξ[3]/2) 

# Boundary Conditions on the edges of the cube
function g(ξ)

    # Check whether we are on the boundary 
    cond1 = ξ[1] == 0 || ξ[1] == 1
    cond2 = ξ[2] == 0 || ξ[2] == 1
    cond3 = ξ[3] == 0 || ξ[3] == 1
    
    # This is cheating, but we already know the solution so we use it.
    if cond1 || cond2 || cond3
        return u_sol(ξ)
    else
        error("$ξ is not a boundary point.")
    end

end




#---------------------------------------#
#      Quadric RBFs & Derivatives       #
#---------------------------------------#

ϕ(ζ,ξ) = sqrt(sum(abs2,ζ-ξ) + σ^2)
function Δϕ(ζ,ξ)
    r2 = sum(abs2,ζ-ξ)
    return (2 * r2 + 3*σ^2) / ( r2 + σ^2 )^(1.5)
end


#---------------------------------------#
#           nIterator Points            #
#---------------------------------------#

iter_one_dim = collect(0:1/Δ:1)
iter_one_inn = collect(1/Δ:1/Δ:(1-1/Δ))

# Boundary Points 
itrtr = Iterators.product(iter_one_dim, iter_one_dim)

iter_boundary = unique(vcat(
    reshape([[0,it1,it2] for (it1, it2) in itrtr], :),
    reshape([[1,it1,it2] for (it1, it2) in itrtr], :),
    reshape([[it1,0,it2] for (it1, it2) in itrtr], :),
    reshape([[it1,1,it2] for (it1, it2) in itrtr], :),
    reshape([[it1,it2,0] for (it1, it2) in itrtr], :),
    reshape([[it1,it2,1] for (it1, it2) in itrtr], :)
))

# Inner Grid Points
iter_inn = [[it1, it2, it3] for (it1,it2,it3) in 
    Iterators.product(iter_one_inn, iter_one_inn, iter_one_inn)]
    
iter_inn = reshape(iter_inn, :)

# Validation 
length(iter_one_dim)^3 == length(iter_boundary) + length(iter_inn)

#---------------------------------------#
#        Linear System Generator        #
#---------------------------------------#

# Make the jth column of the linear system 
function col_make(j)

    L = length(iter_boundary)
    N = length(iter_inn)

    Aj = zeros(Float64, L+N)

    # Get point corresponding to column j
    ξj = j <= L ? iter_boundary[j] : iter_inn[j - L]
            
    # Fill in elements of Aj with boundary points 
    for l in 1:L
        ζl = iter_boundary[l]
        Aj[l] = ϕ(ζl, ξj)
    end

    # Fill in elements of Aj with inner points 
    for l in L+1:L+N
        ζk = iter_inn[l - L]
        Aj[l] = Δϕ(ζk, ξj)
    end

    return Aj

end

# Make the constant vector 
function constant_make()

    L = length(iter_boundary)
    N = length(iter_inn)

    gh = zeros(Float64, L+N)

    # Fill in elements for boundary 
    for l in 1:L
        ξl = iter_boundary[l]
        gh[l] = g(ξl)
    end

    # Fill in elements for inner points 
    for l in L+1:L+N
        ξl = iter_inn[l - L]
        gh[l] = h(ξl)
    end

    return gh
   
end

#---------------------------------------#
#    Approximate Solution Constructor   #
#---------------------------------------#

function u_approx(ζ, coef)

    L = length(iter_boundary)
    N = length(iter_inn)

    val = 0

    # Fill in elements for boundary
    for l in 1:L
        ξl = iter_boundary[l]
        val += coef[l] * ϕ(ζ,ξl)
    end

    # Fill in elements for interior
    for l in L+1: L+N
        ξl = iter_inn[l - L]
        val += coef[l] * ϕ(ζ,ξl)
    end

    return val 
end

function u_approxv(ζ, coef)

    L = length(iter_boundary)
    N = length(iter_inn)
    lg = length(ζ)
    out = Array{Float64}(undef,lg)
    for i = 1:lg
	    val = 0
	
	    # Fill in elements for boundary
	    for l in 1:L
	        ξl = iter_boundary[l]
	        val += coef[l] * ϕ(ζ[i],ξl)
	    end
	
	    # Fill in elements for interior
	    for l in L+1: L+N
	        ξl = iter_inn[l - L]
	        val += coef[l] * ϕ(ζ[i],ξl)
	    end
	    out[i] = val
    end
    return out
end

#---------------------------------------#
#             Example Usage             #
#---------------------------------------#
#=
L = length(iter_boundary)
N = length(iter_inn)

# Coefficient Matrix 
A = zeros(Float64, L+N, L+N)
for j in 1:(L+N)
    A[:,j] = col_make(j)
end

# Constant Vector 
b = constant_make()

# Scalar Coefficients for Approximate Solution 
u = A \ b 

# Compare Approximate Solution to Oracle 
ζ = rand(3)
err = abs(u_sol(ζ) - u_approx(ζ, u))
println("At point $ζ, the absolute error is $err.")
=#

#-------------------------------------------------#
#          Consider Unsymetric problem            #
#-------------------------------------------------#

# Make the jth column of the linear system 
function col_make_un(j,points,bound_points)

    L = length(iter_boundary)
    N = length(iter_inn)
    m,_ = size(points)
    n,_ = size(bound_points)
    Aj = zeros(Float64, L+N+m+n)

    # Get point corresponding to column j
    ξj = j <= L ? iter_boundary[j] : iter_inn[j - L]
            
    # Fill in elements of Aj with boundary points 
    for l in 1:L
        ζl = iter_boundary[l]
        Aj[l] = ϕ(ζl, ξj)
    end
    for l in 1:n
        ζl =  bound_points[l,:]
        Aj[l+L] = ϕ(ζl, ξj)
    end
    # Fill in elements of Aj with inner points 
    for l in 1:N
        ζk = iter_inn[l]
        Aj[l+N+n] = Δϕ(ζk, ξj)
    end
    for l in 1:m
        ζk = points[l,:]
        Aj[L+N+n+l] = Δϕ(ζk, ξj)
    end
    return Aj

end

# Make the constant vector assuming that you provide additional interior points
function constant_make_un(points,bound_points)

    L = length(iter_boundary)
    N = length(iter_inn)
    m,_ = size(points)
    n,_ = size(bound_points)
    gh = zeros(Float64, L+N+m+n)

    # Fill in elements for boundary 
    for l in 1:L
        ξl = iter_boundary[l]
        gh[l] = g(ξl)
    end
    for l in 1:n
        ξl = bound_points[l,:]
        gh[l] = g(ξl)
    end
    # Fill in elements for inner points 
    for l in 1:N
        ξl = iter_inn[l]
        gh[L+l+n] = h(ξl)
    end
    for l in 1:m
        ξl = points[l,:]
        gh[L+N+n+l] = h(ξl)
    end
    return gh
   
end
#=
npoints = 10
rints = rand(npoints,3)
rbounds = rand(npoints,3)
rbounds[1:div(npoints,3),1] .= 0
rbounds[div(npoints,3)+1:2*div(npoints,3),2] .= 0
rbounds[2*div(npoints,3)+1:end,3] .= 0
L = length(iter_boundary)
N = length(iter_inn)
m,_ = size(rints)
n,_ = size(rbounds)
# Coefficient Matrix 
A = zeros(Float64, L+N+m+n, L+N)
for j in 1:(L+N)
    A[:,j] = col_make_un(j,rints,rbounds)
end

# Constant Vector 
b = constant_make_un(rints,rbounds)

# Scalar Coefficients for Approximate Solution 
u = A \ b 

norm(A*u-b)
=#
function test_points(n, sol)
    accum = 0
    for i = 1:n
        xi = rand(3)
        accum += abs(u_sol(xi) - u_approx(xi, sol))
    end
    return accum /n 
end

