function get_setup()
    # Choose DE time-stepper

    # get correct DE functions
   

    if lowercase(type_SC) == lowercase("paulData")
        SC,dist,lags,N,minSC,W_sum,FC_Array = networksetup(c;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
        lags[lags .> 0.0] = lags[lags .> 0.0] .+ constant_delay
        println("mean delay = ", mean(lags[lags .> 0.0]))
    else
        SC,dist,lags,N,minSC,W_sum = networksetup(;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
        FC_Array = []
    end


    ss,NGp,start_adapt,nP,bP,LR,IC,κS,wS,opts,vP,aP,WHISTMAT,d,nRuns,timer,ONES,non_zero_weights,rEcurrent =
    setup(numThreads,nWindows,tWindows;delays=delays,mode=mode,plasticityOpt=plasticityOpt)

    return SC,dist,lags,N,minSC,W_sum,FC_Array,ss,NGp,start_adapt,nP,bP,LR,IC,κS,wS,opts,vP,aP,WHISTMAT,d,nRuns,timer,ONES,non_zero_weights,rEcurrent

end
