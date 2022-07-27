using DifferentialEquations,Parameters,Plots,LinearAlgebra
include("Params.jl")
include("DEFunction.jl")

const PARAMS = LIFParams()
N = 200
W = ones(N,N)
while size(W[W.>0],1)/(N^2) > 0.5
    coord = rand(1:N),rand(1:N)
    if W[coord[1],coord[2]] > 0.0
        W[coord[1],coord[2]] = 0.0
    end
end
const nP = NetworkParams(W,N)
v0 = rand(N)
p=[0]
tspan=(0.0,200.0)
prob = ODEProblem(LIF,v0,tspan,p)
sol = solve(prob,BS3())

heatmap(sol[:,:])