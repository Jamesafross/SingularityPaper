using DifferentialEquations,CircularArrayBuffers,Plots,Random,LinearAlgebra,MKL,SparseArrays,CUDA,CUDA.CUSPARSE,ThreadPinning
include("DEFunctions.jl")
include("functions.jl")
ThreadPinning.mkl_set_dynamic(0)
pinthreads(:compact)

LinearAlgebra.BLAS.set_num_threads(1)
const v_thresh = 10.f0
const v_rest = -v_thresh
const N = 6000
const σ=0.01
const τ = 1.
const η=randn(N)
density = 0.2
W_cpu = zeros(N,N)
generate_matrix!(W_cpu,N,1.0)
const W = CuArray(Float32.(zeros(N,N)))
copyto!(W,W_cpu)
const adapt_time = 0.1f0
const d = Float32.(zeros(N))
const spikeArray = Array{CircularArrayBuffer}(undef,N)
const R = Float32.(zeros(N))
for i = 1:N
    spikeArray[i] = CircularArrayBuffer{Float32}(1,100)
end
const firingTimes = CircularArrayBuffer{Float32}(1,100)
const α = 0.5f0
const κ = 0.9f0/(Nf0)
v0 = Float32.(rand(N))
p=[adapt_time]
tspan=Float32.((0.0,20.0))
prob = ODEProblem(QIFGPU,v0,tspan,p)
@time sol = solve(prob,BS3(),reltol=1e-3)
heatmap(sol[:,:],c=:jet)