include("Collocation_solve.jl")
#Define parameters for intial compliation run
Δ = 3
σ = 1
up = 1e-4
Sketch = "Achlio_2"
s = 20
blocksize = 10
width = 10
max_it = 1000
stop = stop_pa(.9,1.1,.95,.95,up,"Achlio_2")

@time sol = S_Colloc_Solve(Δ, σ, s, Sketch, blocksize, stop, lambda1 = 1, lambda2 = 10, max_it = 100)
@time solkry =  S_Colloc_BLEND(Δ, σ, s, stop,max_it = 100000)
@time M,b = System(Δ,σ)
norm(M*sol[3]-b)
norm(M*solkry[3]-b)
#test on large problem
Δ = parse(Int64,ARGS[1])
Δ = 7
iter_boundary,iter_inn = Prod_bound(Δ)
L = length(iter_boundary)
N = length(iter_inn)
m = L + N
n = L + N
print(m," ",n,"\n")
stop = stop_pa(.9,1.1,.95,.95,up,"Achlio_2")
s = parse(Int64,ARGS[2])
@time sol = S_Colloc_Solve(Δ, σ, s, Sketch, blocksize, stop, lambda1 = 1, lambda2 = 10, max_it = 1000)
len = length(sol[4])
print(len," ",sol[4][len],"\n")
@time solkry = S_Colloc_BLEND(Δ, σ, s, stop,max_it = 1000)
@time M,b = System(Δ,σ)
print(norm(M*sol[3]-b),"\n")
print(norm(M*solkry[3]-b),"\n")
