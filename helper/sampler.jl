using LinearAlgebra: sqrt
using Random: randn
using LinearAlgebra
using Hadamard
using Random
using Statistics
using StatsBase
using LoopVectorization
# Generate a random Matrix using the probability 1/2 method proposed by Achlioptas
# s is the sample dimension, m is the dimension of the full space
function Ach(s,m;split =2)
		       if s > 1
        S = rand(m,s)
    else 
        S = rand(m)
    end
    if split == 2
        S = ifelse.(S .<= .5, -1., 1.)
    else
        S = ifelse.(S .<= 1/3, -1., S)
        S = sqrt(3) .* ifelse.(S .<= 2/3, 0., 1.)
    end
    return 1 / sqrt(s) * S
end

# Generate a matrix using the Fast Lindenstrauss Transform proposed by Ailon et al.

function SRHT(s,m)
    D = Diagonal(float(sample([-1,1],m)))
    H = float(hadamard(m))
    H = D*H
    S = H[:,sample(1:m,s)]
    return S * 1 / sqrt(s)
end

function FJLT(s,m,q)
           D = Array{Int8}(sample([-1,1],m))
           H = hadamard(m)
           H = float(D .* H)
           P = sqrt(1/q)  * randn(m,s)
           q1 = 1 - q
           P = Array{Float64}(ifelse.(rand(m,s) .<= q1, 0, P))
           P = LinearAlgebra.BLAS.gemm('N','N',1/sqrt(m),H,P)
           return P * 1 / sqrt(s)
       end

function gaus(s,m)
    if s == 1
        return randn(m)
    else 
        return 1 / sqrt(s) * randn(m,s)
    end
end

# Function that will generate the random sketching matrix of a specified size
function random_mat(m,s,rand_samp = "gauss", q = 1/4)
    if rand_samp == "gauss"
        return gaus(s,m)
    elseif rand_samp == "uniform"
        return (2 .*rand(m,s) .- 1)
    elseif rand_samp  == "SRHT"
        return SRHT(s,m)
    elseif rand_samp == "Achlio_2"
        return Ach(s,m,split = 2)
    elseif rand_samp == "Achlio_3"
        return Ach(s,m,split = 3)
    elseif rand_samp == "FJLT"
    	return   FJLT(s,m,q)
    elseif rand_samp == "Uniform"
        return Matrix{Float64}(I, m, m)[:, sample(1:m, s, replace = false, ordered = true)]
    end
end


# 2 option achlio
function Ach2!(S)
    s,m = size(S)
    sc = 1/sqrt(m)
    for j = 1:m 
        for i = 1:s 
            if S[i,j] < .5 
                S[i,j] = -1. * sc
            else
                S[i,j] = 1. * sc 
            end
        end
    end
end

#three option achlio
function Ach3!(S)
    s,m = size(S)
    sc = 1/sqrt(m) * sqrt(3)
    for j = 1:m 
        for i = 1:s 
            if S[i,j] < 1/3 
                S[i,j] = -1. * sc
            elseif S[i,j] >= 1/3 && S[i,j] < 2/3
                S[i,j] = 0
            else
                S[i,j] = 1. * sc
            end
        end
    end
end


function Ach!(S;split = 2)
    s,m = size(S)
    rand!(S)
    if split == 2
        Ach2!(S)
    else
        Ach3!(S)
    end
end

#Assumes that column dimension dictates sampling
function GenSke!(A::Array{Float64,2},rand_samp)
    n,m = size(A)
    s = min(n,m)
    d = max(m,n) 
    if rand_samp == "gauss"
        scal = 1/sqrt(s)
        randn!(A)
        BLAS.scal!(d*s,scal,A,1)
        #elseif rand_samp == "uniform"
         #   return (2 .*rand(m,s) .- 1)
    elseif rand_samp  == "SRHT"
        SRHT(s,d)
    elseif rand_samp == "Achlio_2"
        Ach!(A,split = 2)
    elseif rand_samp == "Achlio_3"
        Ach!(A,split = 3)
        #elseif rand_samp == "FJLT"
         #   return  FJLT(s,m,q)
    elseif  rand_samp == "Identity"
        A = Matrix{Float64}(I,s,n)
    end
end


function SRHT_Sketch!(A,s,C)
    m,n = size(A)
    indx = sort(sample(1:m,s,replace = false))
    signs = bitrand(n) * 2 .- 1
    #cols = bitrand(n)
    A .*= signs
    C = fwht(A,1)[indx,:]
    return C
end


#this is a correct implementation of the fast walsh Hadamard trnaform it takes and vector or a matrix as a input
# INPUT: 
# a- One or two dimensional array
# r- a vector if you wish to scale by some modification of an identity matrix\
# OUTPUT:
# Modifies the array A 
function fwht!(a::Array{Float64};r=[1],s_d = 1)
    ti = 0
    ti += @allocated if typeof(a) == Array{Float64,1}
        h = 1
        ln = length(a)
        if rem(log(2,ln),1) != 0.
            print("Vector must be a power of 2 in dimension\n")
            return
        end
        s2 = sqrt(2)
        while h < ln
            inc = 2*h
            for i = 1:inc:ln
                for j = i:(i+h-1)
                    x = a[j]
                    y = a[j+h]
                    a[j] = x + y
                    a[j+h] = x - y
                end
            end
            h *= 2
        end
        if s_d != 1
            scal = (sqrt(ln/s_d)/2^(log2(ln) / 2))
            BLAS.scal!(ln,scal,a,1)
        else
            scal = 1/2^(log2(m) / 2)
            BLAS.scal!(ln,scal,a,1)
        end
    elseif typeof(a) == Array{Float64,2}
        m,n = size(a)
        if rem(log(2,m),1) != 0.
            print("Matrix must be a power of 2 in row dimension\n")
            return
        end
        if length(r) == m
            @turbo a .*= r
        end
        s2 = sqrt(2)
        h = 1
        while h < m
            inc = 2*h
            for k = 1:n
	            for i = 1:inc:m
	                for j = i:(i+h-1)
	                    x = a[j,k]
	                    y = a[j+h,k]
	                    a[j,k] = x + y
	                    a[j+h,k] = x - y
	                end
	            end
            end
            h *= 2
        end
        if s_d != 1
            scal = (sqrt(m/s_d)/2^(log2(m) / 2))
            BLAS.scal!(n*m,scal,a,1)
        else
            scal = 1/2^(log2(m) / 2)
            BLAS.scal!(n*m,scal,a,1)
        end
    else
        print("Need to input either a float matrix or float vector\n")
        return
    end
end

function rbits(n)
    return float.(bitrand(n) .* 2 .- 1)
end

function hadamard_sketch!(A,S; type_s = "row")
    m,n = size(A)
    ms,ns = size(S)
    if type_s == "row"
        if ns != m
            print("Need the row dimension of A to match the col dimension of S")
            return
        end
        low_pow = rem(log(2,m),1)
        if low_pow != 0.
            pow = 2^Int64(div(log(2,m),1) + 1)
            B = zeros(pow,n)
            Av = view(A,1:m,1:n)
            B[1:m,1:n] = Av
            signs = rbits(pow)
            fwht!(B,r=signs,s_d =ms)
            indx = sample(1:pow,ms, replace = false)
            S .= view(B,indx,:)
        else
            signs = rbits(m)
            fwht!(A,r=signs,s_d=ms)
            indx = sample(1:m,ms, replace = false)
            S .= view(A,indx,:)
        end 
    #=else
        if ms != n
            print("Need the row dimension of A to match the row dimension of S")
            return
        end 
        low_pow = rem(log(2,n),1)
        if low_pow != 0.
            pow = 2^Int64(div(log(2,n),1))
            B = zeros(pow,)
            Av = view(A,1:m,1:n)
            B[1:m,1:n] = Av
            signs = bitrand(pow) * 2 .- 1
            fwht!(B,r=signs)
            indx = sample(1:pow,ms, replace = false)
            S = view(B,indx,:)
        else
            signs = bitrand(m) * 2 .- 1
            fwht!(A,r=signs)
            indx = sample(1:pow,ms, replace = false)
            S = view(A,indx,:)
        end=#
    end

end

#Hadamard sampling assuming that you have already preallocated blocks of the appropiate size

function randflipandcpy!(A,B,bits)
    m,n = size(A)
    mb,nb = size(B)
    @inbounds for i = 1:n
        @turbo for j = 1:m
            B[j,i] = A[j,i] * bits[j]
        end
    end
    @inbounds for i = 1:nb
        @turbo for j = (m+1):mb
            B[j,i] = 0
        end
    end
end

function randflipandcpy_vec!(a,bv,bits)
    m = length(a)
    n = length(bv)
    @turbo for i = 1:m
        bv[i] = a[i] * bits[i]
    end
    @turbo for i = m+1:n
        bv[i] = 0
    end
end


function BHSK!(A,S,B,d,indx,start_i, end_i)
    ma,_ = size(A)
    ms,_ = size(S)
    #mb,nb = size(B)
    randflipandcpy!(A,B,d)
    fwht!(B,s_d = ms)
    S[:,start_i:end_i] .= view(B,indx,start_i:end_i)
end 

#Doing sampling assuming padding matrix B has already been created
function BHSK_p!(S,B,d,indx,start_i, end_i)
    ms,_ = size(S)
    fwht!(B,r= d, s_d = ms)
    S[:,start_i:end_i] .= view(B,indx,:)
end 

function BHSK_vec!(b,bs,bsv,d,indx)
    ms = length(bs)
    randflipandcpy_vec!(b,bsv,d) 
    fwht!(bsv,s_d = ms)
    bs .= view(bsv,indx)
end
