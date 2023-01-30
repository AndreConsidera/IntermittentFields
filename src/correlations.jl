function expkernel(r::Union{AbstractArray{<:Real},Real}, ξ::Real)
    @.(exp(-r^ξ))
end
function piecewisekernel(r::Union{AbstractArray{<:Real},Real}, ξ::Real)
    @.((1. - r^ξ) * (r<1)) 
end
function logkernel(r::Union{AbstractArray{<:Real},Real})
    @.(log(1/r) * (r < 1))
end
