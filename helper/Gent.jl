using LoopVectorization
#--------------------------------------------#
#           Decompostion functions           #
#--------------------------------------------#

#function to update the R matrix in a function neeeded for gentleman's
function updateRx!(x,R,ncol,i,xi,c,s)
    for k = (i+1):ncol
        tmp = x[k]
        x[k] -= xi * R[i,k]
        R[i,k] = c * R[i,k] + s * tmp
    end
end

# Implementation of Modified Version of Gentleman's that utilizes only the R matrix. The biggest change is that it does not perform the diivision in the algorithm that would make all diagonal entries one.
# INPUT: 
#   R- Upper triangular amtrix
#   x - row of matrix being decomposed
#   D - Vector that stores non-sqrt scaling factor for R
#   ncols - number of columns in matrix being decomposed
#OUTPUT: Modifies the matrix R
function gentleman!(R ::Array{Float64,2}, x::Array{Float64,1}, D::Array{Float64,1}, ncols::Int64; weight = 1.)
    w = weight
    @inbounds for i = 1:ncols
        if w == 0
            return
        end
            xi = x[i]
        if x[i] != 0 
            tmp = D[i]
            wxi = w * xi
            dp = tmp + wxi * xi
            c = tmp / dp
            s = wxi / dp
            w = c * w
            D[i] = dp
            if i != ncols
                updateRx!(x,R,ncols,i,xi,c,s) 
            end
            
            R[i,i] = 1 
        end
    end
end

#function that copies a column from a matrix to a vector
function copyVecFromMat!(nrows,col,into,from)
    @turbo for i = 1:nrows
        into[i] = from[i,col]
    end
end

#Gentleman applied to a block of the matrix column-wise
# INPUT:
#   A: Matrix block used to update R and D
#   R: Upper triangular matrix being updated
#   D: Vector to be updated
# OUTPUT: Modifies R and D
function BQR!(A, R::Array{Float64,2}, D::Array{Float64,1}, x::Array{Float64,1})
    m,n = size(A)
    for i = 1:n
        copyVecFromMat!(m,i,x,A)
        gentleman!(R,x,D,m)
    end
end

#----------------------------------------------#
#           R Solving functions                #
#----------------------------------------------#


# Transpose upper triangle solve for the gentleman output Needs both the R Matrix and D vector
# INPUT: 
#   R- Upper triangular matrix that will be transposed
#   D- Scaling vector from gentleman'scaling
#   b - constant vector of system
#   x - solution vector
# Output: modifies solution vector
function TraUpperTriSol!(R::Array{Float64,2},D::Array{Float64,1},b::Array{Float64,1},x::Array{Float64,1})
    m,_ = size(R)
    @inbounds for j=1:m
        accum = b[j]
        @inbounds for i = 1:(j-1)
            accum -=  x[i] * R[i,j]
        end
        x[j] = accum 
    end
    @turbo for j = 1:m
        x[j] /= D[j]
    end
end

#Solves a general upper triangular system
# INPUT: 
#   R- Upper triangular matrix that will be transposed
#   D- Scaling vector from gentleman'scaling
#   b - constant vector of system
#   x - solution vector
# Output: modifies solution vector
function UpperTriSol!(R::Array{Float64,2},b::Array{Float64,1},x::Array{Float64,1})
    _,m = size(R)
    @inbounds for j=m:-1:1
        accum = b[j]
        @inbounds for i = m:-1:(j+1)
            accum -=  x[i] * R[j,i]
        end
        x[j] = accum
    end
end

#Function that generates update for Gower Richtarik based using constant vector and matrix R
#INPUT: 
#   R: Upper triangular matrix from gentleman's
#   D: Vector from gentleman's
#   b: Constant vector for the solution
#   up: Update vector
# OUTPUT: modifies the update vector

function GowUp!(R::Array{Float64,2}, D::Array{Float64,1}, b::Array{Float64,1}, inter::Array{Float64}, up::Array{Float64,1})
    ti = 0
    ti += @allocated TraUpperTriSol!(R,D,b,inter)
    ti += @allocated UpperTriSol!(R,inter,up)
    return ti
end

