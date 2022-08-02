include("./functions/NextGen_InitSetup.jl")



println("Base Number of Threads: ",numThreads," | BLAS number of Threads: ", BLASThreads,".")
nWindows = 40
tWindows = 50
type_SC = "pauldata" #sizes -> [18, 64,140,246,503,673]
size_SC = 140
densitySC= 0.3
delay_digits=5
plasticity="on"
mode="rest+perturb"  #(rest,stim or rest+stim)
n_Runs=1
nFC = 1
nSC = 1
eta_0E = -14.7
kappa = 0.0504
delays = "on"
multi_thread = "on"
constant_delay = 0.002
if numThreads == 1
    global multi_thread = "off"
end

NGp = NextGen2PopParams2(η_0E = eta_0E,κ=kappa)

c = 12000

plotdata = "true"

if lowercase(type_SC) == lowercase("paulData")
    const SC,dist,lags,N,minSC,W_sum,FC,missingROIs = networksetup(c;digits=delay_digits,nSC=nSC,nFC=nFC,type_SC=type_SC,N=size_SC,density=densitySC)
    lags[lags .> 0.0] = lags[lags .> 0.0] .+ constant_delay
    println("mean delay = ", mean(lags[SC .> 0.0]))
else
    const SC,dist,lags,N,minSC,W_sum = networksetup(;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
    FC_Array = []
end


const solverStruct =
setup(numThreads,nWindows,tWindows;delays=delays,plasticity=plasticity,NGp=NGp)

if multi_thread == "off"
    nextgen(du,u,h,p,t) = nextgen_unthreads(du,u,h,p,t)
else
    nextgen(du,u,h,p,t) = nextgen_threads(du,u,h,p,t)
end

run_nextgen()

time_per_second = solverStruct.timer.meanIntegrationTime/(2*tWindows)
print(time_per_second)








