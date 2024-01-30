using Statistics: BLAS
using LinearAlgebra
using DataFrames
using CSV
using MatrixDepot
using Random
include("./solver.jl")
include("./sampler.jl")

#Function designed to record all iteration values 
function SketchAndSolve(A::Matrix, b::Vector, xsol::Vector, Binv::Matrix, s::Int64,
    Sketch::String; max_it = 10000, width = 15, FileName = "qewqds", anorm = -1, bnorm = -1)
    m,n = size(A)
    x = zeros(m)
    #Holds the true residuals
    r_true_n = zeros(max_it)
    #Holds the true average residual over same window of sketched moving average
    av_true_res_n = zeros(max_it)
    # Holds moving average residual
    mov_res_n = zeros(max_it)
    # Holds observed residual
    res_ob_n = zeros(max_it)
    # Holds error term
    e_n = zeros(max_it) 
    for i = 1:max_it
        S = random_mat(m, s, Sketch)
        t_res = A * x - b
        o_res = gower_Richtarik!(x,A,b,S,Binv)
        err = xsol - x
        res_ob_n[i] = (BLAS.dot(o_res,o_res))
        r_true_n[i] = (BLAS.dot(t_res,t_res))
        e_n[i] = BLAS.dot(err,err)
        mov_res_n[i] = move_av(res_ob_n,i,width)
        av_true_res_n[i] = move_av(r_true_n,i,width)
    end
    if anorm < 0
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n)
    else
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n, 
            A_Norm = fill(anorm,max_it), B_Norm = fill(bnorm,max_it))
    end 
    if FileName != "qewqds"
        CSV.write(FileName,df)
    end
    return df
end


function SketchAndSolveVM(A::Matrix, b::Vector, xsol::Vector, Binv::Matrix, s::Int64,
    Sketch::String; max_it = 10000, width = 15, FileName = "qewqds", anorm = -1, bnorm = -1, xstart = 0)
    m,n = size(A)
    if typeof(xstart) <: Number 
        x = zeros(m)
    else 
        x = deepcopy(xstart)
    end
    #Holds the true residuals
    r_true_n = zeros(max_it)
    #Holds the true average residual over same window of sketched moving average
    av_true_res_n = zeros(max_it)
    # Holds moving average residual
    mov_res_n = zeros(max_it)
    # Holds observed residual
    res_ob_n = zeros(max_it)
    # Holds error term
    e_n = zeros(max_it) 
    M = Matrix{Float64}(I, m, m)
    for i = 1:max_it
        S = random_mat(m, s, Sketch)
        #idx = sample([1:10, 11:20, 21:30, 31:40, 41:50, 51:60, 61:70, 71:80, 81:90, 91:100, 101:110, 111:120, 121:130, 131:140, 141:150, 151:160, 161:170, 171:180, 181:190, 191:200, 201:210, 211:216],1)
        #S = M[:,idx[1]]
        t_res = A * x - b
        o_res = gower_Richtarik!(x,A,b,S,Binv)
        #SA = S'A
        #Sr = S't_res
        #xp = Krylov.lsqr(SA, Sr)
        #x = xp[1]
        #o_res = SA * x - S'b
        err = xsol - x
        res_ob_n[i] = (BLAS.dot(o_res,o_res))
        r_true_n[i] = (BLAS.dot(t_res,t_res))
        e_n[i] = (BLAS.dot(err,err))
        mov_res_n[i] = move_av(res_ob_n,i,width)
        av_true_res_n[i] = move_av(r_true_n,i,width)
    end
    if anorm < 0
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n)
    else
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n, 
            A_Norm = fill(anorm,max_it), B_Norm = fill(bnorm,max_it))
    end 
    if FileName != "qewqds"
        CSV.write(FileName,df)
    end
    return x,df
end

function test(direct, matrix, row, column, sample_method, sample_size, width, max_it, seed;i = 0)
    m = row
    n = column
    Random.seed!(3123)
    A = Array(matrixdepot(matrix, m))
    Random.seed!(seed)
    xsol = rand(n)
    b = A * xsol

    if m != n
        B = A'A
    elseif m == n
        B = Matrix{Float64}(I,m,n)
    end
    smaxB = opnorm(B)
    smaxA = opnorm(A)
    Binv = inv(B)
    if i == 0
    	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,".csv")
    else
	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,"_",i,".csv")
    end
	SketchAndSolve(A, b, xsol, Binv, sample_size, sample_method,
        max_it = max_it, width = width, FileName = Filename, anorm = smaxA, bnorm = smaxB);
end



#Function designed to record all iteration values 
function SketchAndSolvedw(A::Matrix, b::Vector, xsol::Vector, Binv::Matrix, s::Int64,
    Sketch::String; max_it = 10000, width = 15, FileName = "qewqds", anorm = -1, bnorm = -1)
    m,n = size(A)
    x = zeros(m)
    #Variable that will indicate when the resiudal norm squared was less than 1
    thres = -1
    #Holds the true residuals
    r_true_n = zeros(max_it)
    #Holds the true average residual over same window of sketched moving average
    av_true_res_n = zeros(max_it)
    # Holds moving average residual
    mov_res_n = zeros(max_it)
    # Holds observed residual
    res_ob_n = zeros(max_it)
    # Holds error term
    e_n = zeros(max_it) 
    for i = 1:max_it
        S = random_mat(m, s, Sketch)
        t_res = A * x - b
	o_res = gower_Richtarik!(x,A,b,S,Binv)
        err = xsol - x
        res_ob_n[i] = BLAS.dot(o_res,o_res)
        r_true_n[i] = BLAS.dot(t_res,t_res)
        e_n[i] = BLAS.dot(err,err)
	if res_ob_n[i] < 1
		if thres == -1
			thres = i 
		end
		width = ifelse(i - thres + 1 < 100, i - thres + 1, 100)
	end
        mov_res_n[i] = move_av(res_ob_n,i,width)
        av_true_res_n[i] = move_av(r_true_n,i,width)
    end
    if anorm < 0
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n)
    else
        df = DataFrame(Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n, 
            A_Norm = fill(anorm,max_it), B_Norm = fill(bnorm,max_it))
    end 
    if FileName != "qewqds"
        CSV.write(FileName,df)
    end
    return df
end
#Test function that calls the sketch and solve function with different widths of moving average depending if the mean of the observed residual was less than 1
function testdw(direct, matrix, row, column, sample_method, sample_size, width, max_it, seed;i = 0)
    m = row
    n = column
    Random.seed!(3123)
    A = Array(matrixdepot(matrix, m))
    Random.seed!(seed)
    xsol = rand(n)
    b = A * xsol

    if m != n
        B = A'A
    elseif m == n
        B = Matrix{Float64}(I,m,n)
    end
    smaxB = opnorm(B)
    smaxA = opnorm(A)
    Binv = inv(B)
    if i == 0
    	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,".csv")
    else
	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,"_",i,".csv")
    end
	SketchAndSolvedw(A, b, xsol, Binv, sample_size, sample_method,
        max_it = max_it, width = width, FileName = Filename, anorm = smaxA, bnorm = smaxB);
end

function SketchAndSolveRep(A::Matrix, b::Vector, xsol::Vector, Binv::Matrix, s::Int64,
    Sketch::String; max_it = 1000, width = 1, FileName = "qewqds", anorm = -1, bnorm = -1,nreps = 100)
    m,n = size(A)
    x = zeros(n)
    xup = zeros(n)
    xold = zeros(n)
    iter = zeros(nreps * max_it * width) 
    rep = zeros(nreps * max_it * width) 
    wid = zeros(nreps * max_it * width)
    #Holds the true residuals
    r_true_n = zeros(nreps * max_it * width)
    #Holds the true average residual over same window of sketched moving average
    av_true_res_n = zeros(nreps * max_it * width)
    # Holds moving average residual
    mov_res_n = zeros(nreps * max_it * width)
    # Holds observed residual
    res_ob_n = zeros(nreps * max_it * width)
    # Holds error term
    e_n = zeros(nreps * max_it * width) 
    for i = 1:max_it
    	for j = 1:nreps
		x  = deepcopy(xold)
		for k = 1:width
        		S = random_mat(m, s, Sketch)
        		t_res = A * x - b
        		o_res = gower_Richtarik!(x,A,b,S,Binv)
        		if j == 1 && k == 1
				xup = deepcopy(x) 
			end
			err = xsol - x
        		res_ob_n[(i-1) * nreps * width + (j  - 1) * width +  k] = sqrt(BLAS.dot(o_res,o_res))
        		r_true_n[(i-1) * nreps * width + (j - 1) * width + k] = sqrt(BLAS.dot(t_res,t_res))
        		e_n[(i-1) * nreps * width + (j - 1) * width + k] = BLAS.dot(err,err)
        		mov_res_n[(i-1) * nreps * width + (j - 1) * width + k] = move_av(res_ob_n,i,width)
        		av_true_res_n[(i-1) * nreps * width + (j - 1) * width + k] = move_av(r_true_n,i,width)
			iter[(i-1) * nreps * width + (j - 1) * width + k] = i 
			rep[(i-1) * nreps * width + (j - 1) * width + k] = j 
			wid[(i-1) * nreps * width + (j - 1) * width + k] = k 
		end
    	end
	xold = deepcopy(xup)
    end
    if anorm < 0
        df = DataFrame(it = iter, reps = rep, ob = wid, Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n)
    else
        df = DataFrame(it = iter, reps = rep, ob = wid, Obs_Res = res_ob_n, True_Res = r_true_n, Error = e_n,  Mov_Obs_Res = mov_res_n, Mov_True_Res = av_true_res_n, 
            A_Norm = fill(anorm,max_it*nreps*width), B_Norm = fill(bnorm,max_it*nreps*width))
    end 
    if FileName != "qewqds"
        CSV.write(FileName,df)
    end
    return df
end

function testrep(direct, matrix, row, column, sample_method, sample_size, width, max_it, seed;i = 0,nreps = 1000)
    m = row
    n = column
    Random.seed!(3123)
    A = Array(matrixdepot(matrix, m))
    Random.seed!(seed)
    xsol = rand(n)
    b = A * xsol

    if m != n
        B = A'A
    elseif m == n
        B = Matrix{Float64}(I,m,n)
    end
    smaxB = opnorm(B)
    smaxA = opnorm(A)
    Binv = inv(B)
    if i == 0
    	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,".csv")
    else
	Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,"_",i,".csv")
    end
	SketchAndSolveRep(A, b, xsol, Binv, sample_size, sample_method,
        			max_it = max_it, width = width, FileName = Filename, anorm = smaxA, bnorm = smaxB);
end

function testls(direct, row, column, sample_method, sample_size, width, max_it, seed)
    m = row
    n = column
    h = zeros(m+n)
    L = Matrix{Float64}(I,m+n,m+n)
    Random.seed!(seed)
    A = rand(m,n)
    xsol = ones(n)

    h[n+1:m+n] = A*xsol + randn(m)
    L[1:m,m+1:m+n] = A
    L[m+1:m+n,1:m] = A'
    L[m+1:m+n,m+1:m+n] = zeros(n,n)

    B = Matrix{Float64}(I,m+n,m+n)
    
    Filename = string(direct,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,".csv")
    y=SketchAndSolve(L, h,  zeros(n+m), B, sample_size, sample_method,
        max_it = max_it, width = width, FileName = Filename);
        
    return y
end

function testbig(direct, matrix, row, column, sample_method, sample_size, width, max_it, seed)
    m = row
    n = column
    Random.seed!(seed)
    A = matrixdepot(matrix, m)
    xsol = rand(n)
    b = A * xsol

    Binv = Matrix{Float64}(I,m,n)

    Filename = string(direct,matrix,"_",row,"_",column,"_",sample_method,"_",sample_size,"_",width,".csv")
    s=SketchAndSolve(A, b, xsol, Binv, sample_size, sample_method,
        max_it = max_it, width = width, FileName = Filename);
    return s
end
