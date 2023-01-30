abstract type AbstractNoise end

struct UnitaryWhiteNoise<:AbstractNoise
    wk::Vector{ComplexF64}
    n::Integer

    function UnitaryWhiteNoise(n::Integer)
        wk = (1/2^0.5) * (rand(Normal(0,1), n) + im * rand(Normal(0,1), n))
        
        new(wk, n)
    end

end

function UnitaryWhiteNoise(n::Integer, performdiv::Bool) 
    if performdiv == true
        UnitaryWhiteNoise(div(n, 2) + 1) 
    else 
        UnitaryWhiteNoise(n)
    end
end

struct GmcNoise<:AbstractNoise
    wr::Vector{Real}
    wk::Vector{ComplexF64}

    function GmcNoise(ker::SingularKernel, g1::UnitaryWhiteNoise, g2::UnitaryWhiteNoise, alpha::Real)
        g1_x = irfft(g1.wk, length(ker.r))
        σ_sq = 2 * sum(abs.(ker.Lk[2:end]).^2)
        wr = g1_x .* exp.(alpha .* realization(ker, g2) .- alpha^2 .* σ_sq)
        wk = rfft(wr)

        new(wr, wk)
    end
    
    function GmcNoise(ker::SingularKernel, g1::UnitaryWhiteNoise, g2::UnitaryWhiteNoise, alpha::Real, Pforward::FFTW.rFFTWPlan, Pinv::AbstractFFTs.ScaledPlan)
        g1_x = Pinv * g1.wk
        σ_sq = 2 * sum(abs.(ker.Lk[2:end]).^2)
        wr = g1_x .* exp.(alpha .* realization(ker, g2, Pinv) .- alpha^2 .* σ_sq)
        wk = Pforward * wr
        
        new(wr, wk)
    end

end
