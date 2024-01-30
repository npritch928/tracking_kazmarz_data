include("sampler.jl")

set.seed!(123)

x = rand(512)
vals = zeros(1000,3)
k = 0

for p = 1:3
    global  k += 1
      for i = 1:1000
          S = random_mat(512,p,"FJLT")
          vals[i,k] = abs(norm(S'x)^2 - norm(x)^2)
      end
end


  nx2 = norm(x)^2
  exce = zeros(1001)
  i = 0
  for delt = 0:.01:10
    global i += 1
    exce[i] = sum(vals[:,3] .> delt * nx2) / 1000
  end
  
a = 0:.01:10

nx2 = norm(x)^2
exce = zeros(1001)
i = 0
for delt = 0:.01:10
  global i += 1
  exce[i] = sum(vals[:,2] .> delt * nx2) / 1000
end

a = 0:.1:10