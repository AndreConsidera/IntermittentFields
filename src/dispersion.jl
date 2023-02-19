function dispersion(Np::Integer, dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, p0::AbstractVector{Float64}, sizefunc::Function )
    dx = 2 * pi/N    
    eta = 2 * dx
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
    exittime = zeros(Np);
    exittime2 = zeros(Np);
    exittime4 = zeros(Np);
    exittime8 = zeros(Np);
    exittime16 = zeros(Np);
    exittime32 = zeros(Np);
    exittime64 = zeros(Np);
    finalp = zeros(length(p0), Np);
    finalp2 = zeros(length(p0), Np);
    finalp4 = zeros(length(p0), Np);
    finalp8 = zeros(length(p0), Np);
    finalp16 = zeros(length(p0), Np);
    finalp32 = zeros(length(p0), Np);
    finalp64 = zeros(length(p0), Np);
    for k in 1:Np
        p = p0
        t = 0
        # other times
        t2 = tmax
        t4 = tmax
        t8 = tmax
        t16 = tmax
        t32 = tmax
        t64 = tmax
        # other p's
        p2 = nothing 
        p4 = nothing 
        p8 = nothing 
        p16 = nothing
        p32 = nothing
        p64 = nothing
        for i in 1:maxiter
            #L2 norm with cross correlation
            R = sizefunc(p)

            if R >= lam/2
                t2 = min(t2, i * dt)
                if isnothing(p2) 
                    p2 = p     
                end
            end
            if R >= lam/4
                t4 = min(t4, i * dt) 
                if isnothing(p4) 
                    p4 = p     
                end
            end
            if R >= lam/8
                t8 = min(t8, i * dt)
                if isnothing(p8) 
                    p8 = p     
                end
            end
            if R >= lam/16
                t16 = min(t16, i * dt)
                if isnothing(p16) 
                    p16 = p     
                end
            end
            if R >= lam/32
                t32 = min(t32, i * dt)
                if isnothing(p32) 
                    p32 = p     
                end
            end
            if R >= lam/64
                t64 = min(t64, i * dt)
                if isnothing(p64) 
                    p64 = p     
                end
            end
            t = t + dt
            if R >= lam break end

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p = p + uinterp(p) * dt^0.5 + (2 * kappa * dt)^0.5 * rand(Normal(0,1), length(p))
        end
        exittime[k] = t
        exittime2[k] = t2
        exittime4[k] = t4
        exittime8[k] = t8
        exittime16[k] = t16
        exittime32[k] = t32
        exittime64[k] = t64
        
        finalp[:, k] = p
        finalp2[:, k] = p2
        finalp4[:, k] = p4
        finalp8[:, k] = p8
        finalp16[:, k] = p16
        finalp32[:, k] = p32
        finalp64[:, k] = p64
        
    end
    return exittime, exittime2, exittime4, exittime8, exittime16, exittime32, exittime64, finalp, finalp2, finalp4, finalp8, finalp16, finalp32, finalp64 
end
