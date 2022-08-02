include("./functions/NextGen_InitSetup.jl")

mutable struct times_save
    time
    network_size
    network_density
 end

size_SC_vec = collect(40:20:200)
densitySC_vec = [0.2]

times_save_array = Array{times_save}(undef,size(size_SC_vec,1),size(densitySC_vec,1))

counter_i = 1
numThreads = Threads.nthreads()
if numThreads > 1
    LinearAlgebra.BLAS.set_num_threads(1)
end
BLASThreads = LinearAlgebra.BLAS.get_num_threads()
pinthreads(:compact)


println("Base Number of Threads: ",numThreads," | BLAS number of Threads: ", BLASThreads,".")


for i in size_SC_vec
    counter_j = 1
    for j in densitySC_vec
        global nWindows = 4
        global tWindows = 30
        global type_SC = "generated"
        global size_SC = i
        global densitySC= j
        global delay_digits=3
        global plasticityOpt="off"
        global mode="rest"
        global n_Runs=1
        global eta_0E = -14.19
        global kappa = 0.505
        global delays = "on"
        global multi_thread = "on"
        include("$HOMEDIR/COMP_paper/NextGeneration/RunNextGenBase.jl")

        times_save_array[counter_i,counter_j] = times_save(time_per_second,size_SC,densitySC)
        counter_j += 1

    end
    global counter_i += 1
end

times_plot = zeros(size(size_SC_vec,1))

for i = 1:size(size_SC_vec,1)
    times_plot[i] = times_save_array[i,1].time
end

save("$HOMEDIR/COMP_paper/NextGeneration/times_save.jld","times_save_array",times_save_array)


