module IntermittentFields
    using FFTW
    using Random
    using Statistics
    using Distributions
    using Interpolations

    export normtorus, normevans
    export expkernel, piecewisekernel, logkernel
    export CovarianceCorrelation, CovarianceKernel, SingularCorrelation, SingularKernel
    export UnitaryWhiteNoise, GmcNoise
    export realization, doexittime2, quickrealization, doexittime3, dispersion
    export structurefunc, meanenergy, δv 
    export DiscreteCallback, DiscreteCallbackSet, VectorCallback, VectorCallbackSet

    include("norms.jl")
    include("correlations.jl")
    include("kernels.jl")
    include("callbacks.jl")
    include("noise.jl")
    include("helpers.jl")
    include("stats.jl")
    include("dispersion.jl")
end # end of module