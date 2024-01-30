using LinearAlgebra
using LoopVectorization
using Hadamard
using IterativeSolvers
using StatsBase
include("./Collocation_solve.jl")

ncoords = 2
Δ = 4


desc = rand(ncoords,2)
function GenFixed(Δ)
    iter_one_dim = collect(0:100/Δ:100)
    iter_one_inn = collect(100/Δ:100/Δ:(100-100/Δ))
    N = length(iter_one_inn)^3
    L = (Δ+1)^3 - N
    # Boundary Points 
    itrtr = Iterators.product(iter_one_dim, iter_one_dim)

    iter_bound = unique(vcat(
        reshape([[0,it1,it2] for (it1, it2) in itrtr], :),
        reshape([[100,it1,it2] for (it1, it2) in itrtr], :),
        reshape([[it1,0,it2] for (it1, it2) in itrtr], :),
        reshape([[it1,100,it2] for (it1, it2) in itrtr], :),
        reshape([[it1,it2,0] for (it1, it2) in itrtr], :),
        reshape([[it1,it2,100] for (it1, it2) in itrtr], :)
    ))

    # Inner Grid Points
    iter_inn = [[it1, it2, it3] for (it1,it2,it3) in 
        Iterators.product(iter_one_inn, iter_one_inn, iter_one_inn)]
    
   
    return (Δ+1)^3, vcat(reshape(iter_inn, :), iter_bound)
end

function Genpoints!(ncoords,bound,coord,desc,lb = 1/6, le = 1/6)
    # bound = Array{Bool}(undef,ncoords)
    # coord = rand(ncoords,3)
    # desc = rand(ncoords,2)
    @inbounds for i = 1:ncoords
        bound[i] = 0
        if desc[i,1] < lb
            bound[i] = 1
            #Consider placement on one of the 6 sides
            if desc[2] < 1/6
                coord[i,1] = 0
            elseif  desc[i,2] < 1/3
                coord[i,2] = 0
            elseif desc[i,2] < .5 
                coord[i,3] = 0
            elseif desc[2] <  2/3
                coord[i,1] = 1
            elseif  desc[i,2] < 5/6
                coord[i,2] = 1
            else
                coord[i,3] = 1
            end
        elseif desc[i,1] < lb + le
            bound[i] = 1
            #Consider placement on one of the 12 edges
            if desc[2] < 1/12
                coord[i,1] = 0
                coord[i,2] = 0
            elseif  desc[i,2] < 2/12
                coord[i,1] = 0
                coord[i,3] = 0
            elseif desc[i,2] < 3/12 
                coord[i,3] = 0
                coord[i,2] = 0
            elseif desc[2] < 4/12
                coord[i,1] = 1
                coord[i,2] = 0
            elseif desc[2] <  5/12
                coord[i,1] = 0
                coord[i,2] = 1
            elseif  desc[i,2] < 6/12
                coord[i,1] = 1
                coord[i,3] = 0
            elseif  desc[i,2] < 7/12
                coord[i,1] = 0
                coord[i,3] = 1
            elseif desc[i,2] < 8/12
                coord[i,3] = 1
                coord[i,2] = 0
            elseif desc[i,2] < 9/12
                coord[i,3] = 0
                coord[i,2] = 1
            elseif desc[2] < 10/12
                coord[i,1] = 1
                coord[i,2] = 1
            elseif  desc[i,2] < 11/12
                coord[i,2] = 1
                coord[i,2] = 1
            else 
                coord[i,3] = 1
                coord[i,2] = 1
            end
        end
    end
   #return bound,coord
end


function constant_make_ran!(ncoords, bound, coords, b)
    @inbounds for i = 1:ncoords
        ξi = coords[i,:]
        if bound[i] 
            b[i] = g(ξi)
        else
            b[i] = h(ξi)
        end
    end
end

function row_make_ran!(j, fixed_points, coords, bound, ncoords, A; σ = 1)
    ξj = fixed_points[j]
    @inbounds for i = 1:ncoords
        ζl = coords[i,:]
        if bound[i]
            A[j,i] = ϕ(ζl, ξj, σ)
        else
            A[j,i] = Δϕ(ζl, ξj, σ)
        end
    end
end 

function col_make_ran!(j, fixed_points, coords, bound, ncoords, A; σ = 1)
    ξj = fixed_points[j]
    @inbounds for i = 1:ncoords
        ζl = coords[i,:]
        if bound[i]
            A[i,j] = ϕ(ζl, ξj, σ)
        else
            A[i,j] = Δϕ(ζl, ξj, σ)
        end
    end
end 

function GenBlockRan!(nfixed,ncoords,fixed_points, coords, bound, A)
    for j = 1:nfixed
        col_make_ran!(j, fixed_points, coords, bound, ncoords, A)
    end
end

# In this function you need to have copied b to res
function Update_sol!(x,A,R,D,stor,up,res)
    BLAS.gemv!('N',1.,A,x,-1.,res)
    BQR!(A,R,D,stor)
    GowUp!(R,D,res,stor,up)
    BLAS.gemv!('T',-1.,A,up,1.,x)
    return BLAS.dot(res,res)
end
function sample_cord!(coordinates, decision_c, fixed_points, decision_f)
    n = size(fixed_points, 1)
    s = size(coordinates, 1)
    idx = sample(1:n, s)
    coordintes .= fixed_points[idx, :]
    decision_c .= decision_f[idx,:]
end
function Colloc_Solve(nfixed, ncoords, stop; lambda1 = 1, lambda2 = 100, maxit = 1000)
    estsr = Array{Float64,1}(undef,maxit)
    # Holds observed residual
    estsl = Array{Float64,1}(undef,maxit)
    lambda = lambda1
    kp = Inf
    kup = Inf
    rho = zeros(lambda2)
    iota = zeros(lambda2)
    st1,st2 = stopping(s,[lambda1,lambda2],stop)
    #Preallocate the necessary arrays
    #First the fixed point stuff
    fixed_points = Array{Float64,2}(undef,nfixed,3)
    descision_f = Array{Float64,2}(undef,nfixed,2)
    bound_f = Array{Bool}(undef,nfixed)
    # Second the coordinate vectors that will be reused
    coordinates = Array{Float64,2}(undef,ncoords,3)
    descision_c = Array{Float64,2}(undef,ncoords,2)
    bound_c = Array{Bool}(undef,ncoords)
    # Preallocate the tilde A matrix and the tilde b vector 
    A = Array{Float64,2}(undef,ncoords,nfixed) #note here we save in the transpose direction to perform QR so in this case A' b
    b = Array{Float64,1}(undef,ncoords)

    #Generate the storage for the QR
    R = zeros(ncoords,ncoords)
    D = zeros(ncoords)
    stor = Array{Float64,1}(undef,ncoords)
    res = Array{Float64,1}(undef,ncoords)
    BLAS.blascopy!(ncoords,b,1,res,1)
    up = zeros(ncoords) 
    x = zeros(nfixed)
    rand!(fixed_points)
    rand!(descision_f)
    #Allocate the fixed points
    Genpoints!(nfixed,bound_f,fixed_points,descision_f)
    pos = 1
    for i = 1:maxit
        # Fill in the random points 
        #rand!(coordinates)
        #rand!(descision_c)
        sample_cord!(coordinates, decision_c, fixed_points, decision_f)
        #Generate first set of random points
        Genpoints!(ncoords,bound_c,coordinates,descision_c)

        #Fill the tilde A and tilde b
        constant_make_ran!(ncoords, bound_c, coordinates, b)

        GenBlockRan!(nfixed,ncoords,fixed_points,coordinates,bound_c,A)
        res = BLAS.blascopy!(ncoords,b,1,res,1)
        nrm2 = Update_sol!(x,A,R,D,stor,up,res)
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
        rhok = w_mean(rho,lambda)
        iotak = w_mean(iota,lambda)
        estsr[i] = rhok
        estsl[i] = iotak
 
        BLAS.scal!(ncoords * ncoords,0.,R,1)
        BLAS.scal!(ncoords,0.,D,1)
        if rhok < stop.up && iotak < stop_p
            return m,n,x,estsr[1:(i-1)], estsl[1:(i-1)]
        end
    end
    return m,n,x,estsr, estsl
end


function TRhoEst(A,b,x,fixed_points,nfixed,coordinates,descision_c,bound_c;nsamples = 200)
    rho = 0
    for i = 1:nsamples
       # rand!(coordinates)
       # rand!(descision_c)
        sample_cord!(coordinates, decision_c, fixed_points, decision_f)
        #Generate first set of random points
        Genpoints!(ncoords,bound_c,coordinates,descision_c)

        #Fill the tilde A and tilde b
        constant_make_ran!(ncoords, bound_c, coordinates, b)

        GenBlockRan!(nfixed,ncoords,fixed_points,coordinates,bound_c,A)
        BLAS.gemv!('N',1.,A,x,-1.,b)
        rho += BLAS.dot(b,b)
    end
    return rho/nsamples
end



function Colloc_Solve_rho_est(Δ, ncoords, stop; lambda1 = 1, lambda2 = 100, maxit = 1000)
    be = time()
    nfixed = (Δ+1)^3
    estsr = Array{Float64,1}(undef,maxit)
    # Holds observed residual
    estsl = Array{Float64,1}(undef,maxit)
    estst = Array{Float64,1}(undef,maxit)
    width = Array{Float64,1}(undef,maxit) 
    lambda = lambda1
    kp = Inf
    kup = Inf
    rho = zeros(lambda2)
    trho = zeros(lambda2)
    iota = zeros(lambda2)
    st1,st2 = stopping(s,[lambda1,lambda2],stop)
    #Preallocate the necessary arrays
    #First the fixed point stuff
    #fixed_points = Array{Float64,2}(undef,nfixed,3)
    #descision_f = Array{Float64,2}(undef,nfixed,2)
    #bound_f = Array{Bool}(undef,nfixed)
    # Second the coordinate vectors that will be reused
    coordinates = Array{Float64,2}(undef,ncoords,3)
    descision_c = Array{Float64,2}(undef,ncoords,2)
    bound_c = Array{Bool}(undef,ncoords)
    # Preallocate the tilde A matrix and the tilde b vector 
    A = Array{Float64,2}(undef,ncoords,nfixed) #note here we save in the transpose direction to perform QR so in this case A' b
    b = Array{Float64,1}(undef,ncoords)

    #Generate the storage for the QR
    R = zeros(ncoords,ncoords)
    D = zeros(ncoords)
    stor = Array{Float64,1}(undef,ncoords)
    res = Array{Float64,1}(undef,ncoords)
    BLAS.blascopy!(ncoords,b,1,res,1)
    up = zeros(ncoords) 
    x = zeros(nfixed)
    #rand!(fixed_points)
    #rand!(descision_f)
    #Allocate the fixed points
    nfixed,fixed_points = GenFixed(Δ)
    pos = 1
    for i = 1:maxit
        # Fill in the random points 
        #rand!(coordinates)
        #rand!(descision_c)
        sample_cord!(coordinates, decision_c, fixed_points, decision_f)
        #Generate first set of random points
        Genpoints!(ncoords,bound_c,coordinates,descision_c)
        
        tnrm2 = TRhoEst(A,b,x,fixed_points,nfixed,coordinates,descision_c,bound_c)
        #Fill the tilde A and tilde b
        constant_make_ran!(ncoords, bound_c, coordinates, b)

        GenBlockRan!(nfixed,ncoords,fixed_points,coordinates,bound_c,A)
        res = BLAS.blascopy!(ncoords,b,1,res,1)
        nrm2 = Update_sol!(x,A,R,D,stor,up,res)
        
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
        trho[pos] = tnrm2
        iota[pos] = nrm4
        #update rho and iota
        trhok = w_mean(trho,lambda) 
        rhok = w_mean(rho,lambda)
        iotak = w_mean(iota,lambda)
        estsr[i] = rhok
        estsl[i] = iotak
        estst[i] = trhok
        width[i] = lambda
 
        BLAS.scal!(ncoords * ncoords,0.,R,1)
        BLAS.scal!(ncoords,0.,D,1)
        te = time()
        if rhok < stop.up && iotak < stop_p || te - be > 432000
            return nfixed,ncoords,x,estsr[1:(i-1)], estsl[1:(i-1)],estst[1:(i-1)],width[1:(i-1)]
        end
    end
    return nfixed,ncoords,x,estsr, estsl,estst,width
end


function MatrixExpectationSV(ncoords,Δ,ntrials)
    nfixed = (Δ+1)^3
    #Preallocate the necessary arrays
    #First the fixed point stuff
    #fixed_points = Array{Float64,2}(undef,nfixed,3)
    #descision_f = Array{Float64,2}(undef,nfixed,2)
    #bound_f = Array{Bool}(undef,nfixed)
    # Second the coordinate vectors that will be reused
    coordinates = Array{Float64,2}(undef,ncoords,3)
    descision_c = Array{Float64,2}(undef,ncoords,2)
    bound_c = Array{Bool}(undef,ncoords)
    # Preallocate the tilde A matrix and the tilde b vector 
    A = Array{Float64,2}(undef,ncoords,nfixed)
    B = zeros(nfixed,nfixed)

    #rand!(fixed_points)
    #and!(descision_f)
    #Genpoints!(nfixed,bound_f,fixed_points,descision_f)
    nfixed,fixed_points = GenFixed(Δ)
    max = 0
    for i = 1:ntrials
        rand!(coordinates)
        rand!(descision_c)
        Genpoints!(ncoords,bound_c,coordinates,descision_c)
        GenBlockRan!(nfixed,ncoords,fixed_points,coordinates,bound_c,A)
        maxP = opnorm(A)#maximum(A)
        max = max > maxP ? max : maxP
        B .+= A'A
    end
    B ./= ntrials
    S = svd(B).S
    return S[1], S[end], max
end
