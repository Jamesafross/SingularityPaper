function setup(numThreads,nWindows,tWindows;delays="on",plasticity="on",n_Runs=1,NGp=NextGen2PopParams2())

    W=zeros(N,N)
    W.=SC


    stimOpt = "off"
    stimNodes = [39]
    Tstim = [60,90]
    stimStr = -5.
    stimWindow = 20
    adapt = "off"
    synapses = "1stOrder"
    tPrev = 0.01
    timeAdapt = 0.02
    startWindow = 5
    perturbWindow = 10
    mode = "rest"

   
    
    if lowercase(plasticity) == "on"
        nSave = Int((nWindows-(startWindow-1))*10*tWindows) + 2
    else 
        nSave = 2
    end

    
    κSEEv = ones(N)*NGp.κSEE
    κSIEv = ones(N)*NGp.κSIE
    κSEIv = ones(N)*NGp.κSEI
    κSIIv = ones(N)*NGp.κSII
    κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]

    init0 = round.(make_init_conds(NGp,N),digits=4)

    nP = NetworkParameters(W, dist,lags, N, minSC,W_sum)
    bP = balloonModelParameters()

    LR = 0.0001 # learning rate adaptation
    IC = init(init0)
    κS = weights(κSEEv, κSIEv, κSEIv, κSIIv, κSUM )
    wS = weightSave(zeros(N,nSave),zeros(N,nSave),zeros(N,nSave),zeros(N,nSave),zeros(N,N,nSave),zeros(N,N,nSave),1)
   
    WHISTMAT = zeros(N,N)
    d = zeros(N)
    nRuns = n_Runs
    timer = TimerStruct(0.,0.,0.)
    stimOpts = StimOptions(stimOpt,stimWindow,stimNodes,stimStr,Tstim)
    runOpts = RunOptions(tWindows,nWindows,perturbWindow,mode)
    solverOpts =  SolverOptions(delays,plasticity,adapt,synapses)
    runPars = RunParameters(tPrev,timeAdapt,0,WHISTMAT,d)
    adaptPars = adaptParameters(LR,startWindow,10.01,IC.u0[1:N])

    solverStruct = NextGenSolverStruct(NGp,nP,bP,IC,κS,wS,stimOpts,runOpts,solverOpts,runPars,adaptPars,nRuns,timer)

   
  

 


    return solverStruct

end

