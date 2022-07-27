function run_nextgen()
    @unpack NGp,nP,bP,IC,κS,wS,stimOpts,runOpts,solverOpts,runPars,adaptPars,nRuns,timer = solverStruct
    IC.u0 = make_init_conds(NGp,nP.N)  + 0.1*rand(8*nP.N)

    BOLD_BOTH = NextGenSol_both(0.0,0.0)
    

    for run_mode in runOpts.modeSwitch


        nP.W .= SC
        solverStruct.κS.κSEEv = ones(nP.N)*NGp.κSEE
        solverStruct.κS.κSIEv = ones(nP.N)*NGp.κSIE
        solverStruct.κS.κSEIv = ones(nP.N)*NGp.κSEI
        solverStruct.κS.κSIIv = ones(nP.N)*NGp.κSII
        solverStruct.κS.κSUM  = κS.κSEEv[1]+κS.κSIEv[1]+κS.κSEIv[1]+κS.κSIIv[1]
        solverStruct.IC.u0 = make_init_conds(NGp,N)  + 0.1*rand(8N)
        
        
        solverStruct.wS.κSEEv[:,1] = solverStruct.κS.κSEEv
        solverStruct.wS.κSIEv[:,1] = solverStruct.κS.κSIEv
        solverStruct.wS.κSEIv[:,1] = solverStruct.κS.κSEIv
        solverStruct.wS.κSIIv[:,1] = solverStruct.κS.κSIIv
        solverStruct.wS.count = 2
        
        solverStruct.stimOpts.stimOpt = setstim

        println("Running model ... ")
        @time out = nextgen_model_windows_run()

        BOLD_OUT=[]
        for ii = 1:runOpts.nWindows
                if ii == 1
                    BOLD_OUT= out[:,:,ii]
                else
                    BOLD_OUT = cat(BOLD_OUT,out[:,:,ii],dims=2)
                end
        end


        if run_mode == normal
            global BOLD_BOTH.BOLD_normal = BOLD_OUT
        else
            global BOLD_BOTH.BOLD_perturbed = BOLD_OUT
        end

        
    
        




    end
    if lowercase(mode) == "rest"
        global BOLD = BOLD_REST
    elseif lowercase(mode) == "stim"
        global BOLD = BOLD_STIM
    elseif lowercase(mode) == "rest+stim"
        global BOLD = NextGenSol_both(BOLD_REST.BOLD_rest,BOLD_STIM.BOLD_stim)
    end

end