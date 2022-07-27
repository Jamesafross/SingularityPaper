function Θ(t)
    if t < 0
        return 0.0
    else
        return 1.0
    end
end

function s(t,α)
    return (α^2)*t*exp(-α*t)*Θ(t)
end

function g(t,α,κ,spikeTimes)
    sumg = 0.0
    for i = 1:length(spikeTimes)
        sumg += s(t-spikeTimes[i],α)
    end
    sumg = sumg*κ
    return sumg
end


function generate_matrix!(A,N,density)
    num = Int(round(N*density))

    A[:,1:num] .= 1.

    for i = 1:N
        shuffle!(@view A[i,:])
    end
    return nothing
end

function init_connectivity_mat!(A,N,density)
    generate_matrix(A,N,density)

end

function shuffle_matrix!(A,N)
    @inbounds for i = 1:N
        shuffle!(@view A[i,:])
    end
    return nothing
end



function gpu_thresh!(y,v_thresh,v_rest)
    index = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride = gridDim().y * blockDim().y
    @inbounds for i = index:stride:length(y)
        if y[i] > v_thresh
            y[i] = v_rest
        end  
    end
    return
end

function get_thresholded!(y, v_thresh,v_rest)
    numblocks = ceil(Int, length(y)/256)
    CUDA.@sync begin
        @cuda threads=256 blocks=numblocks gpu_thresh!(y, v_thresh,v_rest)
    end
end

function get_R!(R,spikeArray,v,v_thresh,v_rest,t)
    @inbounds for i = 1:N
        if v[i] > v_thresh
            v[i] = v_rest
            push!(spikeArray[i],t)
            R[i] = g(t,α,κ,spikeArray[i])
        end
    end
end

function get_dv!(dv,v,d)
    @inbounds Threads.@threads for i = 1:N
        dv[i] = τ*(v[i]^2 + 100. + d[i])
    end
end

    
