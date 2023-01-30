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
    export realization, doexittime, quickrealization, doexittime3
    export structurefunc, meanenergy, δv 

    include("norms.jl")
    include("correlations.jl")
    include("kernels.jl")
    include("noise.jl")
    include("helpers.jl")
    include("stats.jl")
end # end of module