@testset "Pikal Microcollapse Tests" begin
	# Test ParamObjMC constructor and indexing
	params_tuple = (
		(nothing, 1.0, 0.1, 1.0),
		(x -> 1.0, 1.0, 1.0),
		(t -> 1.0, t -> 1.0),
		(1.0, 1.0, 1.0, x -> 1.0, x -> 1.0, 1.0)
	)
	param_obj = ParamObjMC(params_tuple)
	@test param_obj[1] == (nothing, 1.0, 0.1, 1.0)
	@test param_obj[2] == (x -> 1.0, 1.0, 1.0)
	@test param_obj[3] == (t -> 1.0, t -> 1.0)
	@test param_obj[4] == (1.0, 1.0, 1.0, x -> 1.0, x -> 1.0, 1.0)
    @test_throws BoundsError(param_obj[0])
    @test_throws BoundsError(param_obj[5])

	# Test lyo_mc_dae! function
	du = zeros(4)
	u = [1.0, 273.15, 0.0, 1.0]
	t = 0.0
	lyo_mc_dae!(du, u, param_obj, t)
	@test !any(isnan.(du))

	# Test ODEFunction creation
	prob = ODEProblem(lyo_mc_dae_f, u, (0.0, 1.0), param_obj)
	sol = solve(prob, Tsit5())
	@test sol.retcode == :Success
end