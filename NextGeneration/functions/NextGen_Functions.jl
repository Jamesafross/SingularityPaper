function adapt_global_coupling(hparams,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64},minSC::Float64,W_sum::Vector{Float64})
    @inbounds for ii = 1:N

        @inbounds for jj = 1:N
         
            if W[ii,jj]  > 0.0
                if lags[ii,jj] == 0.0
                    
                     W[ii,jj] += 0.0000001*u[ii]*(u[jj] - h(hparams,t-1.0;idxs=jj))
                else
                     W[ii,jj] += 0.0000001*h(hparams,t-lags[ii,jj];idxs=ii)*(u[jj] - h(hparams,t-1.0;idxs=jj))
                end
                if W[jj,ii] < minSC
                    W[jj,ii] = minSC
                end
            end
        
        end
       
        
        
        W[W .< 0.0] .= 0.0
        W = W./maximum(W)

    end

     @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end

    return W

end

function adapt_local_func(h,hparams,t,κS,NGp,rE,rI,i,N,c;type = "lim")
        @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
        κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
        @unpack κSEEv,κSIEv,κSEIv,κSIIv,κSUM = κS
        
        κSEEv[i] = (κSEEv[i] + c*rE*(rE - h(hparams,t-1.0;idxs = i)))
        κSIEv[i] = (κSIEv[i] + c*rE*(rI - h(hparams,t-1.0;idxs = i+N)))
        κSEIv[i] = (κSEIv[i] + c*rI*(rE - h(hparams,t-1.0;idxs = i)))
        κSIIv[i] = (κSIIv[i] + c*rI*(rI - h(hparams,t-1.0;idxs = i+N)))
        
        limEE = 0.2
        limEI = 0.2
        limIE = 0.2
        limII = 0.2

        if type == "lim"
            if κSEEv[i]  > κSEE + limEE
                κSEEv[i] = κSEE + limEE
            elseif κSEEv[i]  < κSEE - limEE
                κSEEv[i] = κSEE - limEE
            end

            if κSEIv[i]  > κSEI + limEI
                κSEIv[i] = κSEI + limEI
            elseif κSEIv[i]  < κSEI - limEI
                κSEIv[i] = κSEI - limEI
            end

            if κSIEv[i]  > κSIE + limIE
                κSIEv[i] = κSIE + limIE
            elseif κSIEv[i]  < κSIE - limIE
                κSIEv[i] = κSIE - limIE
            end

            if κSIIv[i]  > κSII + limII
                κSIIv[i] = κSII + limII
            elseif κSIIv[i]  < κSII - limII
                κSIIv[i] = κSII - limII
            end
        elseif type == "normalised"
            κSEEv[i], κSIEv[i],κSEIv[i], κSIIv[i] = κSUM*[κSEEv[i], κSIEv[i], κSEIv[i], κSIIv[i]]/(κSEEv[i] + κSIEv[i] + κSEIv[i] + κSIIv[i])
        end

    return κSEEv[i],κSIEv[i],κSEIv[i],κSIIv[i]
end


function adapt_local_func2(h,hparams,t,κS,NGp,rE,rI,i,N,c;type = "lim")
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    @unpack κSEEv,κSIEv,κSEIv,κSIIv,κSUM = κS
    
    ΔκSEE = c*rE*(rE - h(hparams,t-1.0;idxs = i))
    ΔκSIE = c*rE*(rI - h(hparams,t-1.0;idxs = i+N))
    ΔκSEI = c*rI*(rE - h(hparams,t-1.0;idxs = i))
    ΔκSII = c*rI*(rI - h(hparams,t-1.0;idxs = i+N))

    κSEEv[i] = (κSEEv[i] + ΔκSEE*adapt_contain_function(ΔκSEE,κSEE))
    κSIEv[i] = (κSIEv[i] + ΔκSIE*adapt_contain_function(ΔκSIE,κSIE))
    κSEIv[i] = (κSEIv[i] + ΔκSEI*adapt_contain_function(ΔκSEI,κSEI))
    κSIIv[i] = (κSIIv[i] + ΔκSII*adapt_contain_function(ΔκSII,κSII))
    
  

return κSEEv[i],κSIEv[i],κSEIv[i],κSIIv[i]
end

function adapt_contain_function(Δκ,κ)
    return (exp(-abs( κ + (Δκ))))
end



function adapt_local_func_nodelay(κS,NGp,rE,rI,i,c;type = "lim")
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    @unpack κSEEv,κSIEv,κSEIv,κSIIv,κSUM = κS
    
    κSEEv[i] = κSEEv[i] + c*rE*(rE)
    κSIEv[i] = κSIEv[i] + c*rE*(rI)
    κSEIv[i] = κSEIv[i] + c*rI*(rE)
    κSIIv[i] = κSIIv[i] + c*rI*(rI)
    
    limEE = 0.2
    limEI = 0.2
    limIE = 0.2
    limII = 0.2

    if type == "lim"
        if κSEEv[i]  > κSEE + limEE
            κSEEv[i] = κSEE + limEE
        elseif κSEEv[i]  < κSEE - limEE
            κSEEv[i] = κSEE - limEE
        end

        if κSEIv[i]  > κSEI + limEI
            κSEIv[i] = κSEI + limEI
        elseif κSEIv[i]  < κSEI - limEI
            κSEIv[i] = κSEI - limEI
        end

        if κSIEv[i]  > κSIE + limIE
            κSIEv[i] = κSIE + limIE
        elseif κSIEv[i]  < κSIE - limIE
            κSIEv[i] = κSIE - limIE
        end

        if κSIIv[i]  > κSII + limII
            κSIIv[i] = κSII + limII
        elseif κSIIv[i]  < κSII - limII
            κSIIv[i] = κSII - limII
        end
    elseif type == "normalised"
        κSEEv[i], κSIEv[i],κSEIv[i], κSIIv[i] = κSUM*[κSEEv[i], κSIEv[i], κSEIv[i], κSIIv[i]]/(κSEEv[i] + κSIEv[i] + κSEIv[i] + κSIIv[i])
    end

return κSEEv[i],κSIEv[i],κSEIv[i],κSIIv[i]
end




function make_init_conds(NGp,N)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    
    params = κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    rE0, rI0, vE0, vI0, gEE0, gEI0, gIE0, gII0 = init_conds_ss(params)
    perturb = 0.0*rand(8*N)

    u0 = zeros(8*N)
    u0[1:N] .= rE0
    u0[N+1:2N] .= rI0
    u0[2N+1:3N] .= vE0
    u0[3N+1:4N] .= vI0
    u0[4N+1:5N] .= gEE0
    u0[5N+1:6N] .= gIE0
    u0[6N+1:7N] .= gEI0
    u0[7N+1:8N] .= gII0

    u0 = u0 + perturb
end



function find_best_fit(R,FC_Array)
    fitRvec = zeros(size(FC_Array,3))
    for i = 1:size(FC_Array,3);
        fitRvec[i] = fitR(R,FC_Array[:,:,i]);
    end

    bestElementsTest = findfirst(x->x==sort(fitRvec, rev=true)[1],fitRvec)
    bestElements = findfirst(x->x==sort(fitRvec, rev=true)[1],fitRvec)

    sort(fitRvec, rev=true)
    for i = 1:size(FC_Array,3)
        currentFit = fitR(R,mean(FC_Array[:,:,bestElements],dims=3)[:,:,1])
        best = findfirst(x->x==sort(fitRvec, rev=true)[i],fitRvec)

        if fitR(R,mean(FC_Array[:,:,cat(bestElements,best,dims=1)],dims=3)[:,:,1]) > currentFit
            bestElements = cat(bestElements,findfirst(x->x==sort(fitRvec, rev=true)[i],fitRvec),dims=1)
        end
    end

    fit = fitR(R,mean(FC_Array[:,:,bestElements],dims=3)[:,:,1])
    return fit,bestElements
    
end



function find_non_zero_weights(W)
    return findall(W .> 0.0)
end

function drdt(rA,vA,gAA,gAB,τA,κVAA,κVAB,ΔA)
    return (1. /τA)*(-gAA*rA - gAB*rA - κVAA*rA - κVAB*rA +2. * rA * vA + (ΔA / (τA*pi)))
end

function dvdt(rA,rB,vA,vB,gAA,gAB,τA,κVAB,VsynAA,VsynAB,η_A0)
    return (1. /τA)*(gAA*(VsynAA - vA) + gAB*(VsynAB - vB) + κVAB*(vB - vA) - (τA^2)*(pi^2) * (rA^2.) +  vA^2. + η_A0) 
end


function R(Z)
    return (1/(τ*pi))*real( (1-conj(Z)) / (1+conj(Z)) )
end

function V(Z)
    return (1/(τ*pi))*real( (1-conj(Z)) / (1+conj(Z)) )
end

function perturbW!(W,a)
    for i = 1:size(W,1)
        for j = 1:size(W,1)
            if W[i,j] > 0.0
                Wij = W[i,j] + a*randn()
                while Wij < 0.0
                    Wij = W[i,j] + a*randn()
                end
                W[i,j] = Wij
            end
        end
    end
end





