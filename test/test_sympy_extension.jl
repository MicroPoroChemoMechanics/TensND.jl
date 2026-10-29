# SymPyPythonCall is a weak dependency: `using TensND` loads neither it nor
# PythonCall, and what calls Python lives in `TensNDSymPyPythonCallExt`. With it
# loaded, the constructors called without coordinates build the same SymPy
# symbols as before (their names and assumptions); without it, they say what to
# load, and a cartesian frame of an unsupported dimension says why it has no
# default.

@testsection "SymPyPythonCall as an extension" begin
    ext = Base.get_extension(TensND, :TensNDSymPyPythonCallExt)
    @test ext isa Module

    # The default coordinates, against the symbols written out.
    x, y, z = symbols("x y z", real = true)
    r, θ, ϕ = symbols("r", positive = true), symbols("θ", real = true), symbols("ϕ", real = true)
    @test getcoords(coorsys_cartesian()) == (x, y, z)
    @test getcoords(coorsys_cartesian(Val(Sym))) == (x, y, z)
    @test getcoords(coorsys_polar()) == (r, θ)
    @test getcoords(coorsys_polar(Val(Sym); canonical = true)) == (r, θ)
    @test getcoords(coorsys_cylindrical()) == (r, θ, symbols("z", real = true))
    @test getcoords(coorsys_cylindrical(Val(Sym))) == getcoords(coorsys_cylindrical())
    @test getcoords(coorsys_spherical()) == (θ, ϕ, r)
    @test getcoords(coorsys_spherical(Val(Sym))) == (θ, ϕ, r)
    spheroidal = (symbols("ϕ", real = true), symbols("p", real = true), symbols("q", positive = true))
    @test getcoords(coorsys_spheroidal()) == spheroidal
    @test getcoords(coorsys_spheroidal(spheroidal)) == spheroidal
    @test getcoords(coorsys_spheroidal(Val(Sym))) == spheroidal
    @test first(init_cartesian()) == (x, y, z)
    @test first(init_cartesian(2)) == symbols("x y", real = true)
    @test first(init_polar()) == symbols("r θ", real = true)
    @test first(init_cylindrical()) == (r, symbols("θ z", real = true)...)
    @test first(init_spherical()) == (symbols("θ ϕ", real = true)..., r)
    @test first(init_rotated()) == symbols("θ ϕ ψ", real = true)
    @test_throws ArgumentError init_cartesian(4)

    # Loaded, the extension adds nothing to a `MethodError`.
    @test sprint(TensND._sympy_defaults_hint, MethodError(coorsys_polar, ()), Any[], Any[]) == ""

    # Without SymPyPythonCall, in a process of its own.
    code = """
    using TensND
    print(any(m -> nameof(m) == :PythonCall, values(Base.loaded_modules)), " ")
    try
        coorsys_spherical()
    catch err
        print(sprint(showerror, err))
    end
    """
    out = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`, String)
    @test startswith(out, "false ")
    @test occursin("using SymPyPythonCall", out)
end
