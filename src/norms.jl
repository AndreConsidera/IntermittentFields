"""
    normevans(x::Real, η::Real)

1D Norm on the torus regularized at scale η.

# Examples
```jldoctest
julia> η = 0.1; normevans(0, η)
0.05

julia> η = 0.; normevans(0, η)
0.0
```

"""
function normevans(x::Real, η::Real)
    xmod = x%(2π)
    xtmp = (xmod < π) * abs(xmod) - (xmod >= π) * abs(2π - xmod)
    tmp = abs(xtmp)
    if η <= 0.
        return tmp
    else
        return 0.5(tmp^2/η + η )*(tmp <= η) + tmp * (tmp > η)
    end 
end