abstract type AbstractKernel end
abstract type AbstractCorrelation end

struct CovarianceCorrelation<:AbstractCorrelation
    c::Function
end 

struct SingularCorrelation<:AbstractCorrelation
    c::Function
end

struct CovarianceKernel<:AbstractKernel
    Lr::Vector{Real}
    Lr_sq::Vector{Real}
    Lk::Vector{ComplexF64}
    cr::Vector{Real}
    ck::Vector{Complex}
    r::Array{Real}
    eta::Real
    rescale::Bool
    cor::CovarianceCorrelation
    ξ::Real
    
    function CovarianceKernel(r::AbstractArray, eta::Real, cor::CovarianceCorrelation, ξ::Real, rescale::Bool)
        reta = normevans(r,eta)
        cr = cor.c(reta, ξ)
        ck = rfft(copy(cr))./length(r)
        if rescale
            # rescale energy
            ck[1] = 0
            ck[end] = 0
            imag(ck[imag(ck).<0]) .= 0
            ck = abs.(ck)/(2. * sum(abs.(ck)))
        end
        Lk = (abs.(copy(ck)).^0.5) 
        Lr_sq = fftshift(irfft(Lk.^2, length(r))) .* length(r)
        Lr = irfft(copy(Lk), length(r)) .* length(r)
        new(Lr, Lr_sq, Lk, cr, ck, r, eta, rescale, cor, ξ)
    end
end

CovarianceKernel(r::AbstractArray, eta::Real, ξ::Real) = CovarianceKernel(r::AbstractArray, eta::Real, CovarianceCorrelation(expkernel), ξ::Real, true)


struct SingularKernel<:AbstractKernel
    Lr::Vector{Real}
    Lr_sq::Vector{Real}
    Lk::Vector{Complex}
    cr::Vector{Real}
    ck::Vector{Complex}
    r::Array{Real}
    eta::Real
    rescale::Bool
    cor::SingularCorrelation
    #ξ::Real
    
    function SingularKernel(r::AbstractArray, eta::Real, cor::SingularCorrelation; rescale::Bool = false)
        reta = normevans(r,eta)
        cr = cor.c(reta)
        ck = rfft(copy(cr))./length(r)
        if rescale
            ck[1] = 0
            ck[end] = 0
            ck[imag(ck).<0].= 0
            ck = @.(abs(ck)/(2*sum(abs(ck))))
        end
        Lk = abs.(copy(ck)).^0.5 
        Lr_sq = fftshift(irfft(Lk.^2, length(r))) .* length(r)
        Lr = irfft(copy(Lk), length(r)) .* length(r)
        new(Lr, Lr_sq, Lk, cr, ck, r, eta, rescale, cor)
    end
end

SingularKernel(r::AbstractArray, eta::Real) = SingularKernel(r::AbstractArray, eta::Real, SingularCorrelation(logkernel))
