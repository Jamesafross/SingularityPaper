using Distributed,SharedArrays
if nprocs() < 20
    addprocs(20)
    println("Number of Workers = ", nworkers())
end

@everywhere begin 
    using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations
    include("./functions/NextGen_InitSetup.jl")
    nWindows = 4
    tWindows = 200
    type_SC = "pauldata"
    size_SC = 140
    densitySC= 0.3
    delay_digits=6
    plasticity="off"
    mode="rest"  #(rest,stim or rest+stim)
    n_Runs=1
    eta_0E = -13
    kappa = 0.2
    delays = "on"
    multi_thread = "on"
    constant_delay = 0.002
    meanFC,missingROI = get_mean_all_functional_data(;ROI=140,type="control")


    if numThreads == 1
        global multi_thread = "off"
    end
    if multi_thread == "off"
        nextgen(du,u,h,p,t) = nextgen_unthreads(du,u,h,p,t)
    else
        nextgen(du,u,h,p,t) = nextgen_threads(du,u,h,p,t)
    end
    c = 7000.
    const SC,dist,lags,N,minSC,W_sum,FC,missingROIs = networksetup(c;digits=delay_digits,type_SC=type_SC,N=size_SC,density=densitySC)
    lags[lags .> 0.0] = lags[lags .> 0.0] .+ constant_delay
        
    const solverStruct =
    setup(numThreads,nWindows,tWindows;delays=delays,mode=mode,plasticity=plasticity)
end








nVec1 = 40
nVec2 = 20
nVec3 = 10
eta_0E_vec = SharedArray(Array(LinRange(-13.5,-15.1,nVec1)))
kappa_vec = SharedArray(Array(LinRange(0.01,0.06,nVec2)))
c_vec = SharedArray(Array(LinRange(5000,14000,nVec3)))

fitArray = SharedArray(zeros(nVec1,nVec2,nVec3))

fitArrayStruct = Array{fitStruct}(undef,nVec1,nVec2,nVec3)



    
@sync @distributed for i = 1:nVec1; 
    for j = 1:nVec2;
        for jj = 1:nVec3
            solverStruct.nP.lags = dist./c_vec[jj]
            solverStruct.nP.lags[solverStruct.nP.lags .> 0.0] = solverStruct.nP.lags[solverStruct.nP.lags .> 0.0] .+ constant_delay
            

            solverStruct.NGp.η_0E = eta_0E_vec[i]
            solverStruct.NGp.κ=kappa_vec[j]
            

            solverStruct.timer.meanIntegrationTime = 0.0 
            run_nextgen()
            
            time_per_second = solverStruct.timer.meanIntegrationTime/tWindows
            print(time_per_second)

         
            if lowercase(type_SC) == "pauldata" 
                modelFC = get_FC(BOLD.BOLD_rest)
                if size(missingROI,1) > 0
                    keepElements = ones(N)
                    for i in missingROI
                        keepElements .= collect(1:N) != i
                    end

                    modelFC = modelFC[keepElements,keepElements]
                end


                FC_fit_to_data_mean = zeros(size(modelFC,3))

                for i = 1:size(modelFC,3)
                    FC_fit_to_data_mean[i] = fit_r(modelFC[:,:,i],meanFC)
                end

            end



            fitArray[i,j,jj] = maximum(FC_fit_to_data_mean)

        end
    end
end


for i = 1:nVec1
    for j = 1:nVec2
        for jj = 1:nVec3
        fitArrayStruct[i,j,jj] = fitStruct(c_vec[jj],eta_0E_vec[i],kappa_vec[j],fitArray[i,j,jj])
        end
    end
end

save("$WORKDIR/NextGen_fitVec.jld","fitArrayStruct",fitArrayStruct)


print(maximum(fitArray))






