using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations
HOMEDIR=homedir()
PROGDIR = "$HOMEDIR/COMP_paper"
WORKDIR="$PROGDIR/NextGeneration"
BALLOONDIR="$PROGDIR/Balloon_Model"
InDATADIR="$HOMEDIR/NetworkModels_Data/StructDistMatrices"
OutDATADIR="$HOMEDIR/NetworkModels_Data/2PopNextGen_Data"
include("$PROGDIR/GlobalFunctions/Global_Headers.jl")
include("$BALLOONDIR/BalloonModel.jl")
include("$WORKDIR/functions/NextGen_Headers.jl")


include("$InDATADIR/getData.jl")

@static if Sys.islinux() 
    using ThreadPinning,MKL
    ThreadPinning.mkl_set_dynamic(0)
    pinthreads(:compact)
end

numThreads = Threads.nthreads()
if numThreads > 1
    LinearAlgebra.BLAS.set_num_threads(1)
end
BLASThreads = LinearAlgebra.BLAS.get_num_threads()