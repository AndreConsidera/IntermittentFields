abstract type AbstractCallback end

struct DiscreteCallback<:AbstractCallback
    condition::Function
    affect::Function
end

struct CallbackSet<:AbstractCallback
    set::Tuple
end