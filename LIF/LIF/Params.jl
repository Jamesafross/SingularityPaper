@with_kw struct LIFParams{R}
    V_r::R = 0.
    V_max::R = 10.
    r::R = 6.
    C::R = 3.
    I::R = 2.
    τ::R = r*C
end

mutable struct NetworkParams
    W
    N
end
