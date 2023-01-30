function normtorus(x::Union{AbstractArray{<:Real},Real})
    xtilde=x.%(2*π )
    tmp=@.((xtilde<π)*abs(xtilde)-(xtilde>=π)*abs(2*π-xtilde))
    return tmp
end

function normevans(x::Union{AbstractArray{<:Real},Real},η::Real)
    tmp = abs.(normtorus(x))
    if η<=0
        return tmp
    else
        return @.((0.5*tmp^2/η +η*0.5 )*(tmp <=η) +tmp *(tmp>η))
    end 
end