using DrWatson, Test
@quickactivate "IntermittentFields"

# Here you include files using `srcdir`
include(srcdir("norms.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "IntermittentFields tests" begin
    @test normevans(0, 0.1) == 0.1/2
    @test normevans(0, 0.1) == normevans(2π, 0.1)
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

