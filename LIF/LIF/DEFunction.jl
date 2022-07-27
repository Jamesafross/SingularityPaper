function LIF(dv,v,p,t)
    @unpack V_r,V_max,r,C,I,τ = PARAMS
    @unpack W,N = nP
    d = zeros(N)
    mul!(d,W,v)
    for i = 1:N
       
        if v[i] > V_max
            v[i] = V_r
        end

       
        dv[i] = (-v[i] + r*I + d[i])/τ 
    end


end