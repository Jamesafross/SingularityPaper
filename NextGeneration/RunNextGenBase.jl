
# Choose DE time-stepper
DDEalg = BS3()
ODEalg = BS3()


# get correct DE functions
if multi_thread == "off"
    nextgen(du,u,h,p,t) = nextgen_unthreads(du,u,h,p,t)
else
    nextgen(du,u,h,p,t) = nextgen_threads(du,u,h,p,t)
end

if lowercase(type_SC) == lowercase("paulData")
    const SC,dist,lags,N,minSC,W_sum,FC_Array = networksetup(c;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
    lags[lags .> 0.0] = lags[lags .> 0.0] .+ constant_delay
    println("mean delay = ", mean(lags[lags .> 0.0]))

else
    const SC,dist,lags,N,minSC,W_sum = networksetup(;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
end

const ss,NGp,start_adapt,nP,bP,LR,IC,κS,wS,opts,vP,aP,WHISTMAT,d,nRuns,timer,ONES,non_zero_weights,rEcurrent =
setup(numThreads,nWindows,tWindows;delays=delays,mode=mode,plasticityOpt=plasticityOpt)

run_nextgen()

time_per_second = timer.meanIntegrationTime/tWindows
