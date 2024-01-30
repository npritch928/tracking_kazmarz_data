using LinearAlgebra
using LoopVectorization
using Hadamard
using IterativeSolvers
include("collocation_make.jl")
include("helper/Gent.jl")
include("helper/sampler.jl")
include("LT_blocking/block_write_read.jl")
include("IterativeSolvers.jl-master/src/IterativeSolvers.jl")


#-------------------------------#
#           Outline             #
#-------------------------------#

# 1. Function that generates a block of the matrix and constant vector
# 2. Sketching Collocation Solver
# 3. Blendenpik Collocation Solver
#---------------------------------------#
#      Quadric RBFs & Derivatives       #
#---------------------------------------#
function abs2_sum_ofdiff(x,y)
    n = length(x)
    if n != length(y)
        print("Vectors not same length\n")
    end
    accum = 0
    @inbounds for i = 1:n
        v = x[i] - y[i]
        accum += v * v
    end
    return accum 
end



function ϕ(ζ,ξ,σ)
    x = abs2_sum_ofdiff(ζ,ξ) + σ^2
    y = sqrt(x)
    return y 
end

function Δϕ(ζ,ξ,σ)
    r2 = abs2_sum_ofdiff(ζ,ξ)
    s2 = σ^2
    return (2 * r2 + 3*s2) / ( r2 + s2 )^(1.5)
end 


#-------------------------------#
#     Matrix Block Generator    #
#-------------------------------#
function col_make!(blj,j,L,N,σ,iter_inn,iter_boundary,A)

    ti = 0
    # Get point corresponding to column j
    ξj = j <= L ? iter_boundary[j] : iter_inn[j - L]
            
    # Fill in elements of Aj with boundary points 
    @inbounds for l in 1:L
        ζl = iter_boundary[l]
        A[l,blj] = ϕ(ζl, ξj,σ)
    end

    # Fill in elements of Aj with inner points 
    @inbounds for l in L+1:L+N
        ζk = iter_inn[l - L]
        A[l,blj] = Δϕ(ζk, ξj,σ)
    end
end

# Make the constant vector 
function constant_make!(iter_boundary, iter_inn, gh)

    L = length(iter_boundary)
    N = length(iter_inn)

    # Fill in elements for boundary 
    @inbounds for l in 1:L
        ξl = iter_boundary[l]
        gh[l] = g(ξl)
    end

    # Fill in elements for inner points 
    @inbounds for l in L+1:L+N
        ξl = iter_inn[l - L]
        gh[l] = h(ξl)
    end

    return gh
   
end

# Function that loads the column values to the block
function GenBlock!(col_start::Int64,col_end::Int64,σ,iter_inn,iter_boundary,B)
    L = length(iter_boundary)
    N = length(iter_inn)
    i = 1
    for j = col_start:col_end
        col_make!(i,j,L,N,σ,iter_inn,iter_boundary,B)
        i+=1
    end
end
#Makes column of a unform sampling function
function col_make!(blj,j,L,N,σ,iter_inn,iter_boundary,A,S)
    s = length(S)
    k = 1
    # Get point corresponding to column j
    ξj = j <= L ? iter_boundary[j] : iter_inn[j - L]
            
    # Fill in elements of Aj with boundary points 
    @inbounds for l in 1:L
        if l == S[k]
            ζl = iter_boundary[l]
            A[k,blj] = ϕ(ζl, ξj,σ)
            k+=1
        end
        if k > s
            break
        end
    end

    # Fill in elements of Aj with inner points 
    @inbounds for l in L+1:L+N
        if l == S[k]
            ζk = iter_inn[l - L]
            A[k,blj] = Δϕ(ζk, ξj,σ)
            k+=1
        end
        if k > s
            break
        end
    end
end
#Makes Constant vector for uniform sampling
function constant_make!(iter_boundary, iter_inn, gh, S)

    L = length(iter_boundary)
    N = length(iter_inn)
    k = 1
    # Fill in elements for boundary 
    @inbounds for l in 1:L
        if l == S[k]
            ξl = iter_boundary[l]
            gh[l] = g(ξl)
            k+=1
        end
    end

    # Fill in elements for inner points 
    @inbounds for l in L+1:L+N
        if l == S[k]
            ξl = iter_inn[l - L]
            gh[l] = h(ξl)
            k+=1
        end
    end

    return gh
   
end

# Function that loads the column values to the block
function GenBlock!(col_start::Int64,col_end::Int64,σ,iter_inn,iter_boundary,B, S)
    L = length(iter_boundary)
    N = length(iter_inn)
    i = 1
    for j = col_start:col_end
        col_make!(i,j,L,N,σ,iter_inn,iter_boundary,B,S)
        i+=1
    end
end
#function that perforrms the subtraction at the end to get the residual
function subres!(res::Array{Float64,1},b::Array{Float64,1})
    m = length(res)
    @turbo for i = 1:m
        res[i] -= b[i]
    end
end
function Prod_bound(Δ)
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
    iter_inn = [[it1, it2, it3] for (it1,it2,it3) in 
    Iterators.product(iter_one_inn, iter_one_inn, iter_one_inn)]
    
    iter_inn = reshape(iter_inn, :)
    return iter_boundary,iter_inn
end
mutable struct stop_pa
    d1::Float32
    d2::Float32
    xi1::Float32
    xi2::Float32
    up::Float64
    Sample_meth::String
end
function stopping(s,lambda,stop::stop_pa)
    if stop.Sample_meth == "gauss" || stop.Sample_meth == "Achlio_2" || stop.Sample_meth == "Achlio_3" 
        C = .23467
        omega = .1127
        eta = 5
    else
        C = .03125 
        omega =  .0625
        eta = 82
    end
    t1 = min((1 - stop.d1)^2 * stop.up,(1 - stop.d1)/omega) / log(2 / stop.xi1)
    t2 = min((stop.d2 - 1)^2 * stop.up,(stop.d2 - 1)/omega) / log(2 / stop.xi2)
    return C * s * lambda * eta * stop.up ./ (1 .+ log.(lambda)) * min(t1,t2)
end

#calcaulates a mean of specified length lambda
function w_mean(u::Array{Float64,1},lambda::Int64)
    if lambda > 5
        accum = 0
        @turbo for i = 1:lambda
            accum += u[i]
        end
    else
        accum = 0
        for i = 1:lambda
            accum += u[i]
        end
    end 
    return accum / lambda
end
#Memory efficient way to apply the matrix to the Solution this can easily be done in parallel
#INPUT:
# S- Random matrix
# B- Block to hold generated matrix
# Bv - view of the block matrix to hold the remaining columns
# up- the update generated from GowUp
# x -  the current solution state
# inter - An array that holds intermediate results will be same size as x
# iter_inn - inner points for boundary equations
# iter_boundary - boundary iterates
# blocksize - self explanatory
# nblocks - the number of blocks to be generated
# cols_rem - the number of left over columns after generating nblocks of size blocksize
# OUTPUT:
# x-Updates the state of x based on all of the above info
function AtSup!(S::Array{Float64,2}, B::Array{Float64,2}, Bv, up::Array{Float64,1},
     x::Array{Float64,1},inter::Array{Float64,1},σ,iter_inn, iter_boundary, blocksize::Int64,nblocks::Int64,cols_rem::Int64)
    col_start = 1
    col_end = 1 + blocksize - 1
    BLAS.gemv!('N',1.,S,up,0.,inter)
    for i = 1:nblocks
        xv = view(x,col_start:col_end)
        GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B) 
        BLAS.gemv!('T',-1.,B,inter,1.,xv)
        col_start = 1 + col_end
        col_end = col_start + blocksize - 1
        #Need to switch postions of x and inter
    end
    #check if number of blocks is odd or not to decide which x or inter goes where
    if cols_rem > 0 
        col_end = col_start + cols_rem - 1
        GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B)
        xv = view(x,col_start:col_end) 
        BLAS.gemv!('T',-1.,Bv,inter,1.,xv)
    end
    return
end
function AtSup!(S::Array{Int64,1}, B::Array{Float64,2}, Bv, up::Array{Float64,1},
     x::Array{Float64,1},inter::Array{Float64,1},σ,iter_inn, iter_boundary, blocksize::Int64,nblocks::Int64,cols_rem::Int64)
    col_start = 1
    col_end = 1 + blocksize - 1
    #inter[S] .= up
    for i = 1:nblocks
        xv = view(x,col_start:col_end)
        GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B,S) 
        BLAS.gemv!('T',-1.,B,up,1.,xv)
        col_start = 1 + col_end
        col_end = col_start + blocksize - 1
        #Need to switch postions of x and inter
    end
    #check if number of blocks is odd or not to decide which x or inter goes where
    if cols_rem > 0 
        col_end = col_start + cols_rem - 1
        GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B,S)
        xv = view(x,col_start:col_end) 
        BLAS.gemv!('T',-1.,Bv,up,1.,xv)
    end
    return
end
#-------------------------------#
#  Sketching Collocation Solver #
#-------------------------------#
Δ = 3
σ = 1
Sketch = "gauss"
s = 30
blocksize = 10
width = 10
max_it = 1000
stop = stop_pa(.9,1.1,.95,.95,1e-10,"gauss")

function S_Colloc_Solve(Δ, σ, s, Sketch, blocksize, stop; lambda1 = 1, lambda2 = 100, max_it = 100000)
    estsr = Array{Float64,1}(undef,max_it)
    # Holds observed residual
    estsl = Array{Float64,1}(undef,max_it)
    # Holds moving average residual
    ti = 0
    si = 0
    oi = 0
    mi = 0
    lambda = lambda1
    kp = Inf
    kup = Inf
    rho = zeros(lambda2)
    iota = zeros(lambda2)
    st1,st2 = stopping(s,[lambda1,lambda2],stop)
    iter_boundary,iter_inn = Prod_bound(Δ)
    L = length(iter_boundary)
    N = length(iter_inn)
    m = L + N
    n = L + N
    x = zeros(m)
    if Sketch == "Uniform" 
    	S = Array{Int64,1}(undef, s)
    else
	S = Array{Float64,2}(undef,m,s)
    end
    AtS = Array{Float64,2}(undef,s,blocksize)
    Stb = Array{Float64,1}(undef,s)
    res = zeros(s)
    Modi = Array{Float64,1}(undef,n)
    b = Array{Float64,1}(undef,m)
    inter = zeros(m)
    B = Array{Float64,2}(undef,m,blocksize)
    stor = zeros(s)
    R = zeros(s,s)
    D = zeros(s)
    constant_make!(iter_boundary, iter_inn, b)
    #b = l
    nblocks = div(N+L,blocksize)
    rem_cols = rem(N+L,blocksize)
    col_start = 1
    col_end = 1 + blocksize - 1
    pos = 1
    #Subsets of matrix required for remaining column calculations
    Bre = view(B,:,1:rem_cols)
    AtSv = view(AtS,:,1:rem_cols)
    ti = 0
    for i = 1:max_it
        if Sketch == "Uniform" 
	        sample!(1:m, S, replace = false, ordered=true) 
            Stb .= b[S] 
	 else	
	    GenSke!(S, Sketch)
	    BLAS.gemv!('T',1.,S,b,0.,Stb)
	 end
        for j = 1:nblocks
            xv = view(x,col_start:col_end)
	        if Sketch == "Uniform"
                GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,AtS,S)
	        else
                GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B)
		        BLAS.gemm!('T','N',1.,S,B,0.,AtS)
	        end
            BQR!(AtS,R,D,stor)
            BLAS.gemv!('N',1.,AtS,xv,1.,res)
            col_start = 1 + col_end
            col_end = col_start + blocksize - 1
        end
        
        #Handle the remained of the columns
        if rem_cols > 0
            col_end = col_start + rem_cols - 1
            xv = view(x,col_start:col_end)
            if Sketch == "Uniform" 
                GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,AtSv,S)
	        else
                GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B)
	    	    BLAS.gemm!('T','N',1.,S,Bre,0.,AtSv)
            end
	        BQR!(AtSv,R,D,stor)
            BLAS.gemv!('N',1.,AtSv,xv,1.,res)
        end    
        subres!(res,Stb)
        nrm2 = BLAS.dot(res,res)
        #use the R matrix to update the solution
        GowUp!(R,D,res,Stb,stor)
        if Sketch == "Uniform"
        # Function to apply A'S to the update (stor) and adds it to x
            AtSup!(S, AtS, AtSv, stor, x,inter,σ,iter_inn, iter_boundary, blocksize, nblocks, rem_cols)
        else
            AtSup!(S, B, Bre, stor, x,inter,σ,iter_inn, iter_boundary, blocksize, nblocks, rem_cols)
        end
        #update Lambda
        if i > 1 && kp == Inf && rho[pos] <= nrm2 
            kp = i
            #The last value of i until we are in the lambda2 regime
            kup = kp + lambda2 - lambda1 - 1
        end
        lambdao = lambda
        if i > kup 
            lambda = lambda2
            stop_p = st2
            pos = rem(i - kp,lambda2) + 1 
        elseif kp == Inf
            lambda = lambda1
            stop_p = st1
            pos = 1
        else 
            lambda = lambda1 + 1 + (i - kp)
            stop_p = stopping(s,lambda,stop)
            pos = lambda 
        end
            
        nrm4 = nrm2^2
        rho[pos] = nrm2
        iota[pos] = nrm4
        #update rho and iota
        rhok = w_mean(rho,lambda)#(rhok * lambdao - rho[pos] + nrm2) / lambda
        iotak = w_mean(iota,lambda)#(iotak * lambdao - iota[pos] + nrm4) / lambda
        estsr[i] = rhok
        estsl[i] = iotak
        #lambdas[i] = pos
        #kps[i] = kp
        BLAS.scal!(s*s,0.,R,1)
        BLAS.scal!(s,0.,D,1)
        BLAS.scal!(s,0.,res,1)
        BLAS.scal!(m,0.,inter,1)
        col_start = 1
        col_end = 1 + blocksize - 1
        #check stopping criterion
        if rhok < stop.up && iotak < stop_p
            #print(i,"\n")
            return m,n,x,estsr[1:(i-1)], estsl[1:(i-1)]#,lambdas,kps#z0 .+ ud, ests
        end
    end
    #A = zeros(Float64, L+N, L+N)
    #GenBlock!(1,L+N,σ,iter_inn,iter_boundary,A)

    return m,n,x,estsr, estsl#,A
end

#function to zero out the lower triangular portion from geqrf
function zero_geqr!(A)
    m,n = size(A)
    @inbounds for i = 1:n
        @turbo for j = (i+1):m
            A[j,i] = 0
        end
    end
end


function S_Colloc_BLEND(Δ, σ, s, stop;max_it = 100000)
    iter_boundary,iter_inn = Prod_bound(Δ)
    L = length(iter_boundary)
    N = length(iter_inn)
    pre = 0
    m = L + N
    n = L + N
    ml2 = log(2,m)
    pow = rem(ml2,1) != 0 ? 2^Int64(div(ml2,1) + 1) : m
    AS = Array{Float64}(undef,s,n)
    A = Array{Float64}(undef,pow,n)
    R = zeros(n,n)
    b = Array{Float64,1}(undef,m)
    inter = Array{Float64,1}(undef,n)
    D = zeros(n)
    tau = Array{Float64}(undef,s)
    Av = view(A,1:m,1:n)
    for i = 1:3
        signs = bitrand(pow) * 2 .- 1
        indx = sample(1:pow,s, replace = false)
        GenBlock!(1,n,σ,iter_inn,iter_boundary,Av)
        #BHSK_p!(AS,A,signs,indx,1,n)
        #LAPACK.geqrf!(AS,tau)
        zero_geqr!(AS)
        #BQR!(AS,R,D,inter)
        #=kap = LAPACK.trcon!('1','U','N',AS)
        if kap > 1.1102230246251565e-15 #the rhs is 5 * eps
            pre = 1
            break
        end=#
    end
    constant_make!(iter_boundary, iter_inn, b)
    #@turbo R .*= sqrt.(D)
    if pre == 1
        #sol = Krylov.lsqr(A,b,N = factorize(AS), ldiv = true)
        sol = IterativeSolvers.lsqr(Av,b,atol = stop.up, maxiter = 100000)#,Pr = factorize(AS))
    else
        sol = IterativeSolvers.lsqr(Av,b,atol = stop.up,maxiter = 100000)
    end 
    return m,n,sol
end


function System(Δ,σ)
    iter_boundary,iter_inn = Prod_bound(Δ)
    L = length(iter_boundary)
    N = length(iter_inn)
    pre = 0
    m = L + N
    n = L + N
    A = Array{Float64,2}(undef,m,n)
    b = Array{Float64,1}(undef,m)
    GenBlock!(1,n,σ,iter_inn,iter_boundary,A)
    constant_make!(iter_boundary,iter_inn,b)
    return A,b
end


function SystemLT(Δ,σ,blocksize,dir,root)
    iter_boundary,iter_inn = Prod_bound(Δ)
    L = length(iter_boundary)
    N = length(iter_inn)
    m = L + N
    n = L + N 
    nblocks = div(n,blocksize)
    cols_rem = rem(n,blocksize)
    A = Array{Float64,2}(undef,m,blocksize)
    Arem = Array{Float64,2}(undef,m,cols_rem)
    for i = 1:nblocks
        start_i = ((i-1) * blocksize + 1)
        end_i = (i * blocksize) 
        GenBlock!(start_i,end_i,σ,iter_inn,iter_boundary,A)
        matWrite_block(dir,root,i,A)
    end
    if cols_rem > 0
        start_i = end_i + 1
        end_i = nblocks * blocksize + cols_rem
        GenBlock!(start_i,end_i,σ,iter_inn,iter_boundary,Arem)
        matWrite_block(dir,root,nblocks + 1,Arem)
    end
end

function Colloc_LSQR_LT!(Δ, σ, dir, memlim;max_it = 10)
    iter_boundary,iter_inn = Prod_bound(Δ)
    m = (Δ+1)^3 
    mem = Int(round(3 * memlim / 4 * 1024 ^ 3))
    colsize = m * 8
    blocksize = min(max(div(mem,colsize),1),m)
    SystemLT(Δ, σ, blocksize, dir, "mat")
    b = Array{Float64,1}(undef,m)
    constant_make!(iter_boundary,iter_inn,b)
    nblocks = div(m,blocksize)
    remcols = rem(m,blocksize)
    x = zeros(m)
    lg = IterativeSolversLT.lsqr!(x,b,dir,"mat",blocksize,nblocks,remcols, maxiter = max_it,log =true)
    return x,lg
end

function Colloc_LSQR_Err!(Δ, σ, b, xsol; max_it = 10)
    iter_boundary,iter_inn = Prod_bound(Δ)
    m = (Δ+1)^3 
    A = Array{Float64,2}(undef,m,m)
    GenBlock!(1,m,σ,iter_inn,iter_boundary,A)
    err = Array{Float64,1}(undef,max_it)
    #constant_make!(iter_boundary,iter_inn,b)
    x = zeros(m)
    for i = 1:max_it
        IterativeSolvers.lsqr!(x,A,b, maxiter = 1)
        err[i] = norm(x-xsol)
    end
    return x,err
end

function S_Colloc_Solve_err(Δ, σ, s, l, xsol, Sketch, blocksize, stop; lambda1 = 1, lambda2 = 100, max_it = 100000)
    estsr = Array{Float64,1}(undef,max_it)
    # Holds observed residual
    estsl = Array{Float64,1}(undef,max_it)
    # Holds moving average residual
    ti = 0
    si = 0
    oi = 0
    mi = 0
    lambda = lambda1
    kp = Inf
    kup = Inf
    rho = zeros(lambda2)
    iota = zeros(lambda2)
    st1,st2 = stopping(s,[lambda1,lambda2],stop)
    iter_boundary,iter_inn = Prod_bound(Δ)
    L = length(iter_boundary)
    N = length(iter_inn)
    m = L + N
    n = L + N
    x = zeros(m)
    S = Array{Float64,2}(undef,m,s) 
    AtS = Array{Float64,2}(undef,s,blocksize)
    Stb = Array{Float64,1}(undef,s)
    res = zeros(s)
    Modi = Array{Float64,1}(undef,n)
    b = Array{Float64,1}(undef,m)
    inter = Array{Float64,1}(undef,m)
    B = Array{Float64,2}(undef,m,blocksize)
    stor = Array{Float64,1}(undef,s)
    R = zeros(s,s)
    D = zeros(s)
    #constant_make!(iter_boundary, iter_inn, b)
    b = l
    nblocks = div(N+L,blocksize)
    rem_cols = rem(N+L,blocksize)
    col_start = 1
    col_end = 1 + blocksize - 1
    pos = 1
    #Subsets of matrix required for remaining column calculations
    Bre = view(B,:,1:rem_cols)
    AtSv = view(AtS,:,1:rem_cols)
    ti = 0
    for i = 1:max_it
        GenSke!(S, Sketch)
        BLAS.gemv!('T',1.,S,b,0.,Stb)
        for j = 1:nblocks
            GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B)
            xv = view(x,col_start:col_end)
            BLAS.gemm!('T','N',1.,S,B,0.,AtS)
            BQR!(AtS,R,D,stor)
            BLAS.gemv!('N',1.,AtS,xv,1.,res)
            col_start = 1 + col_end
            col_end = col_start + blocksize - 1
        end
        
        #Handle the remained of the columns
        if rem_cols > 0
            col_end = col_start + rem_cols - 1
            xv = view(x,col_start:col_end)
            GenBlock!(col_start,col_end,σ,iter_inn,iter_boundary,B)
            BLAS.gemm!('T','N',1.,S,Bre,0.,AtSv)
            BQR!(AtSv,R,D,stor)
            BLAS.gemv!('N',1.,AtSv,xv,1.,res)
        end    
        subres!(res,Stb)
        nrm2 = BLAS.dot(res,res)
        #use the R matrix to update the solution
        GowUp!(R,D,res,Stb,stor)
        # Function to apply A'S to the update (stor) and adds it to x
        AtSup!(S, B, Bre, stor, x,inter,σ,iter_inn, iter_boundary, blocksize, nblocks, rem_cols)
        #update Lambda
        if i > 1 && kp == Inf && rho[pos] <= nrm2 
            kp = i
            #The last value of i until we are in the lambda2 regime
            kup = kp + lambda2 - lambda1 - 1
        end
        lambdao = lambda
        if i > kup 
            lambda = lambda2
            stop_p = st2
            pos = rem(i - kp,lambda2) + 1 
        elseif kp == Inf
            lambda = lambda1
            stop_p = st1
            pos = 1
        else 
            lambda = lambda1 + 1 + (i - kp)
            stop_p = stopping(s,lambda,stop)
            pos = lambda 
        end
            
        nrm4 = nrm2^2
        rho[pos] = nrm2
        iota[pos] = nrm4
        #update rho and iota
        rhok = w_mean(rho,lambda)#(rhok * lambdao - rho[pos] + nrm2) / lambda
        iotak = w_mean(iota,lambda)#(iotak * lambdao - iota[pos] + nrm4) / lambda
        estsr[i] = rhok
        estsl[i] = norm(x - xsol)
        #lambdas[i] = pos
        #kps[i] = kp
        BLAS.scal!(s*s,0.,R,1)
        BLAS.scal!(s,0.,D,1)
        BLAS.scal!(s,0.,res,1)
        col_start = 1
        col_end = 1 + blocksize - 1
        #check stopping criterion
        if rhok < stop.up && iotak < stop_p
            #print(i,"\n")
            return m,n,x,estsr[1:(i-1)], estsl[1:(i-1)]#,lambdas,kps#z0 .+ ud, ests
        end
    end
    #A = zeros(Float64, L+N, L+N)
    #GenBlock!(1,L+N,σ,iter_inn,iter_boundary,A)

    return m,n,x,estsr, estsl#,A
end

function scal_sys!(A,b)
    m = size(A,1)
    for i = 1:m
        n = norm(A[i,:])
        A[i,:] ./= n
        b[i] /= n
    end
end
