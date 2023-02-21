abstract type AbstractCallback end

struct DiscreteCallback<:AbstractCallback
    condition::Function
    affect::Function
end

struct VectorCallback<:AbstractCallback
    condition::Function
    affect::Function
end


struct VectorCallbackSet <:AbstractCallback
    set::Tuple 
end

struct DiscreteCallbackSet <:AbstractCallback
    set::Tuple 
end
