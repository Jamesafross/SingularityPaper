using DifferentialEquations,CircularArrayBuffers,Plots,Random,LinearAlgebra,MKL,SparseArrays,CUDA,CUDA.CUSPARSE
include("DEFunctions.jl")
include("functions.jl")

LinearAlgebra.BLAS.set_num_threads(1)
const v_thresh = 10.
const v_rest = -v_thresh
const N = 2000
const σ=0.01
const W = (init_connectivity_mat(N,0.2))
const d = zeros(N)
const spikeArray = Array{CircularArrayBuffer}(undef,N)
const R = zeros(N)
for i = 1:N
    spikeArray[i] = CircularArrayBuffer{Float64}(1,100)
end
const firingTimes = CircularArrayBuffer{Float64}(1,100)
const α = 0.5
const κ = 0.5/(40)
v0 = ones(N)
p=[0]
tspan=(0.0,200.0)
prob = ODEProblem(QIF,v0,tspan,p)
@time sol = solve(prob,BS3(),reltol=1e-3)
heatmap(sol[:,:])