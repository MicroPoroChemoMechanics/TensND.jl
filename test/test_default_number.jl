# The element type of what is built without one: `Sym` once SymPyPythonCall is
# loaded, `Float64` otherwise. The suite loads SymPyPythonCall; a child process
# that does not load it checks the other case.

@testsection "Default element type" begin
    @test default_number_type() == typeof(one(SymPyPythonCall.Sym)) && isconcretetype(default_number_type())
    @test eltype(tens_J4()) == default_number_type()
    @test tens_K4(2, Float64) == TensISO{2}(0.0, 1.0)
    @test mandel(tens_J4(2, Float64)) ≈ [0.5 0.5 0.0; 0.5 0.5 0.0; 0.0 0.0 0.0]
    @test eltype(levi_civita()) == default_number_type()
    @test CanonicalBasis() isa CanonicalBasis{3, default_number_type()}
    code = """
    using TensND
    T = default_number_type()
    𝕀, 𝕁, 𝕂 = iso_projectors()
    ok = T === Float64 && eltype(tens_Id4()) === Float64 && 𝕁 + 𝕂 == 𝕀 &&
        eltype(𝐞(1)) === Float64 && eltype(levi_civita()) === Float64 && CanonicalBasis() isa CanonicalBasis{3, Float64}
    exit(ok ? 0 : 1)
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`
    @test success(run(ignorestatus(cmd)))
end
