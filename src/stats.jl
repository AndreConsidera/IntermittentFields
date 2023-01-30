function structurefunc(ker::AbstractKernel, m::Integer, p::Integer) 
    r0s = 2 .^[i for i in 0:floor(log2(length(ker.r))) - 1]
    dr=2 * pi/length(ker.r)
    l = r0s * dr
    
    sptmp = zeros(length(r0s), m)
    sp_abstmp = copy(sptmp)
    sp = zeros(length(r0s))
    sp_abs = copy(sp)
    for k in 1:m
        noise = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
        u = realization(ker, noise)
        for (i, h) in enumerate(r0s)
            δv = circshift(u,-h) - u
            sptmp[i, k] = mean(δv.^p)
            sp_abstmp[i, k] = mean(abs.(δv).^p)
        end
    end
    sp = mean(sptmp, dims = 2)[1:end]
    sp_abs = mean(sp_abstmp, dims = 2)[1:end]
    return sp_abs, sp, l
end

function structurefunc(ker::AbstractKernel, m::Integer, p::Integer, α::Real) 
    r0s = 2 .^[i for i in 0:floor(log2(length(ker.r))) - 1]
    dr = 2 * pi/length(ker.r)
    l = r0s * dr
    
    sptmp = zeros(length(r0s), m)
    sp_abstmp = copy(sptmp)
    sp = zeros(length(r0s))
    sp_abs = copy(sp)
    
    logker = SingularKernel(ker.r, ker.eta, SingularCorrelation(logkernel))
    for k in 1:m
        g1 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
        g2 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
        noise = GmcNoise(logker, g1, g2, α)
        u = realization(ker, noise)
        for (i, h) in enumerate(r0s)
            δv = circshift(u,-h) - u
            sptmp[i, k] = mean(δv.^p)
            sp_abstmp[i, k] = mean(abs.(δv).^p)
        end
    end
    sp = mean(sptmp, dims = 2)[1:end]
    sp_abs = mean(sp_abstmp, dims = 2)[1:end]
    return sp_abs, sp, l
end


function meanenergy(ker::AbstractKernel, n::Integer)
    energy = zeros(n)
    m = div(length(ker.r),2) + 1
    for i in eachindex(energy)
        noise = UnitaryWhiteNoise(m)
        u = realization(ker, noise)
    energy[i] = mean(u.^2)
    end
    return  mean(energy)    
end

function meanenergy(ker::AbstractKernel, n::Integer, α::Real)
    energy = zeros(n)
    m = div(length(ker.r),2) + 1
    
    logker = SingularKernel(ker.r, ker.eta, SingularCorrelation(logkernel))
    for i in eachindex(energy)
        g1 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
        g2 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
        noise = GmcNoise(logker, g1, g2, α)
        u = realization(ker, noise)
    energy[i] = mean(u.^2)
    end
    return  mean(energy)    
end

function δv(ker::AbstractKernel, m::Integer)
    r0s = 2 .^[i for i in 0:floor(log2(length(ker.r))) - 1]
    δv = zeros(length(r0s), m, length(ker.r))
    for k in 1:m
        for (i, h) in enumerate(r0s)
            noise = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
            u = realization(ker, noise)
            δvtmp = circshift(u,-h) - u
            δv[i, k, :] = δvtmp
        end
    end
    
    return reshape(δv,(length(r0s), m * length(ker.r)))
end 

function δv(ker::AbstractKernel, m::Integer, α::Real)
    r0s = 2 .^[i for i in 0:floor(log2(length(ker.r))) - 1]
    δv = zeros(length(r0s), m, length(ker.r))
    
    logker = SingularKernel(ker.r, ker.eta)
    for k in 1:m
        for (i, h) in enumerate(r0s)
            g1 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
            g2 = UnitaryWhiteNoise(length(ker.r)÷ 2 + 1)
            noise = GmcNoise(logker, g1, g2, α)
            u = realization(ker, noise)
            δvtmp = circshift(u,-h) - u
            δv[i, k, :] = δvtmp
        end
    end
    
    return reshape(δv,(length(r0s), m * length(ker.r)))
end 
