function QIF(dv,v,p,t)
    println(t)
    copyto!(d,W*R)
    for i = 1:N
        if v[i] > v_thresh
            v[i] = v_rest
            push!(spikeArray[i],t)
            R[i] = g(t,α,κ,spikeArray[i])
        end
    end
    @inbounds Threads.@threads for i = 1:N

        dv[i] = v[i]^2 + 0.1 + d[i]
    end
end

function QIFGPU(dv,v,p,t)
    adapt_time = p[1]
    #println(t)
   # if t > adapt_time
    #    shuffle_matrix!(W_cpu,N)
     #   copyto!(W,W_cpu)
      #  adapt_time += 0.1
    #end
    get_R!(R,spikeArray,v,v_thresh,v_rest,t)
    copyto!(d,W*CuArray(R))
    
    get_dv!(dv,v,d)

     
end

function dw(dv,v,p,t)
    dv[:] .= σ
end
