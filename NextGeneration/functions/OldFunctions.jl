function rij(cij,bsdp)
    return (cij + 1.0)^bsdp
end

function hill(x::Float64,hsdp::Float64,bsdp::Float64)
    return x/(x+hsdp^bsdp) - 1/2
end


function modHeaviside(x)
    if x < 0.0
        return 1.0
    else
        return -1.0
    end
end

function ΔWsdp(cij::Float64,hsdp::Float64,bsdp::Float64)
    return 0.0001*hill(rij(cij,bsdp),hsdp,bsdp)
end

function ΔWgtp(w,dist,cij,agdp,ηgdp)
    return agdp*modHeaviside(w - exp(cij*dist))*ηgdp
end

function adapt_global_coupling_cor(N::Int64,W::Matrix{Float64},dist::Matrix{Float64},minSC::Float64,W_sum::Vector{Float64},HIST::Array{Float64},hsdp,bsdp)
  
@inbounds for ii = 1:N
        
    @inbounds for jj = 1:N 
        if W[jj,ii]  > 0.0
            cij = cor(HIST[ii,:],HIST[jj,:])
            W[jj,ii] += ΔWsdp(cij,hsdp,bsdp)
            if W[jj,ii] < minSC
                W[jj,ii] = minSC
            elseif W[jj,ii] > 0.12
                W[jj,ii] = 0.12
            end
        end
    
    end
   
    if sum(W[:,ii]) != 0.0
    @views W[:,ii] = W_sum[ii].*(W[:,ii]./sum(W[:,ii]))
    end
    W[W .< 0.0] .= 0.0

end

 @inbounds for k=1:N #Maintaining symmetry in the weights between regions
    @inbounds for l = k:N
             W[k,l] = (W[k,l]+W[l,k])/2
             W[l,k] = W[k,l]
    end
end

return W

end
