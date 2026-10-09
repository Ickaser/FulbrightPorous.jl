using DrWatson, Test
@quickactivate "FulbrightPorous"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "FulbrightPorous tests" begin
	include("test_pikal_microcollapse.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
