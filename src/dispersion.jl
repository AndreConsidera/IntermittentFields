function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function )
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end

# with logmessage
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function  )
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
            
            # logmessage
            logmessage(i, maxiter, (length(t[R .< lam])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end
# with λ
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function, λ::Real)
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta

    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= λ break end
            t[R .< λ] = t[R .< λ] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< λ] .= p[:, R .< λ] .+ uinterp.(p[:, R .< λ]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< λ], 2) )
            
            # logmessage
            logmessage(i, maxiter, (length(t[R .< λ])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end

# with DiscreteCallback
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function, cb::DiscreteCallback)
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
            
            #callback
            if cb.condition(p, t, k)
                cb.affect(p, t, k)
            end

            # logmessage
            logmessage(i, maxiter, (length(t[R .< lam])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end

# with DicreteCalbackSet
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function, cb::DiscreteCallbackSet )
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
            
            #DiscreteCallbackSet
            for cb in cb.set
                if cb.condition(p, t, k)
                    cb.affect(p, t, k)
                end
            end

            # logmessage
            logmessage(i, maxiter, (length(t[R .< lam])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end

# with VectorCallback
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function, cb::VectorCallback )
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
            
            #VectorCallback
            for (idx, condition) in enumerate(cb.condition(p, t))
                if condition
                    cb.affect(p, t, idx)
                end
            end

            # logmessage
            logmessage(i, maxiter, (length(t[R .< lam])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end

# with VectorCalbackSet
function dispersion(nindep::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, logmessage::Function, cb::VectorCallbackSet )
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    lam = 1
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = Int(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= lam break end
            t[R .< lam] = t[R .< lam] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< lam] .= p[:, R .< lam] .+ uinterp.(p[:, R .< lam]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< lam], 2) )
            
            #VectorCallbackSet
            for cb in cb.set
                for (idx, condition) in enumerate(cb.condition(p, t))
                    if condition
                        cb.affect(p, t, idx)
                    end
                end
            end

            # logmessage
            logmessage(i, maxiter, (length(t[R .< lam])/ncor)*100)
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime,finalp
end