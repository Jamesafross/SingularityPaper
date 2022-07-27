function generate_network_matrices(N;density=1)
    W=ones(N,N)
    W=W.-diagm(ones(N))
    triSize = (N^2 - N)/2
    numOnes = triSize
    while numOnes/triSize >  density 
        pos_i = rand(1:N)
        pos_j = rand(pos_i:N)
        if W[pos_i,pos_j] == 1
            W[pos_i,pos_j] = 0
            numOnes += -1
        end
    
    end

    #initial shuffle 
    for n = 1:5
        for i = 1:N
            W[i,i+1:N] = shuffle(shuffle(W[i,i+1:N]))
        end
        for j = 2:N
            W[1:j-1,j] = shuffle(shuffle(W[1:j-1,j]))
        end
    end

    while sum(W,dims=2)[sum(W,dims=2) .== 0] != []
        for i = 1:N
            W[i,i+1:N] = shuffle(shuffle(W[i,i+1:N]))
        end
        for j = 2:N
            W[1:j-1,j] = shuffle(shuffle(W[1:j-1,j]))
        end
    end


    lags = W.*(rand(1:100,N,N)/10000)
    lags .= Symmetric(lags)
    W = W.*rand(N,N)
    W = Symmetric(W)
    SC = zeros(N,N)
    SC .= W
    
    for i = 1:N
        SC[:,i] = SC[:,i]/sum(SC[:,i])
    end
    return SC,lags

end



