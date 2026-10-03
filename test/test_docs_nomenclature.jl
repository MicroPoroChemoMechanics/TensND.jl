# The nomenclature of the documentation, docs/nomenclature.toml: the page
# `nomenclature.md` renders it and the hints over the equations read it. Ported
# from MeanFieldHomogenization's test of the same name.
#
# Every page an entry is restricted to exists, a symbol has at most one meaning
# that holds everywhere and at most one on any given page, two entries never
# carry the same name, and the formulas of the pages and of the docstrings follow
# the typography of Notation and conventions.

using TOML

@testsection "Nomenclature of the documentation" begin
    root = pkgdir(TensND)
    docs = joinpath(root, "docs")
    entries = TOML.parsefile(joinpath(docs, "nomenclature.toml"))["symbol"]
    @test length(entries) > 50
    for e in entries
        @test all(k -> haskey(e, k), ("group", "tex", "label", "name"))
        for p in get(e, "pages", String[])
            @test endswith(p, "/") ? isdir(joinpath(docs, "src", p)) : isfile(joinpath(docs, "src", p))
        end
    end

    # One meaning that holds everywhere per form, and one meaning per form on a
    # given page: two scoped entries of the same symbol never share a page.
    form(e) = (e["tex"], get(e, "script", false))
    holds(p, q) = p == q || (endswith(p, "/") && startswith(q, p)) || (endswith(q, "/") && startswith(p, q))
    for f in unique(form.(entries))
        same = [e for e in entries if form(e) == f]
        @test count(e -> isempty(get(e, "pages", String[])), same) <= 1
        scoped = [get(e, "pages", String[]) for e in same if !isempty(get(e, "pages", String[]))]
        for i in eachindex(scoped), j in (i + 1):lastindex(scoped)
            clash = [(p, q) for p in scoped[i] for q in scoped[j] if holds(p, q)]
            isempty(clash) || @error "two meanings of $(f[1]) on the same page" clash
            @test isempty(clash)
        end
    end

    # The label, the name and the unit are shown as HTML over the equations
    # (docs/src/.vitepress/theme/symbol-hints.ts), so they hold no tag but <sub>,
    # <sup>, <b> and <u>, and every one is closed.
    allowed = ("<sub>", "</sub>", "<sup>", "</sup>", "<b>", "</b>", "<u>", "</u>")
    for e in entries, field in ("label", "name", "unit")
        text = get(e, field, "")
        @test all(in(allowed), [m.match for m in eachmatch(r"</?[^>]*>", text)])
        for t in ("sub", "sup", "b", "u")
            @test count("<$t>", text) == count("</$t>", text)
        end
    end

    # One definition per concept: two entries never carry the same name.
    names = [lowercase(strip(e["name"])) for e in entries]
    @test allunique(names)

    # The typography of the formulas, on the pages and in the docstrings, whose
    # formulas the API pages render. What would bring back a second convention,
    # or hide a symbol from the hints, is refused.
    refused = [
        r"\\rm\b" => "{\\rm …}: write \\mathrm{…}",
        r"\\bm\b" => "\\bm: write \\boldsymbol",
        r"\\hat\s*\{\s*\\(underline|mathbf)" => "a hat on a vector",
        r"\\underline\s*\{\s*\\underline" => "a doubly underlined tensor: write \\boldsymbol",
        r"\\mathbf\s*\{?\s*1\s*\}?" => "\\mathbf 1: the identity is \\boldsymbol{1}",
        r"[^\x00-\x7f]" => "a non-ASCII character inside a formula",
    ]
    function formulas(text)
        fences = [m.captures[1] for m in eachmatch(r"```math\n(.*?)```"s, text)]
        rest = replace(text, r"```.*?```"s => " ")
        inline = [m.captures[1] for m in eachmatch(r"``((?:[^`]|`(?!`))+?)``", rest)]
        return vcat(fences, inline)
    end
    sources = Tuple{String, String}[]
    for (r, _, fs) in walkdir(joinpath(docs, "src")), f in fs
        endswith(f, ".md") && !occursin("node_modules", r) && !occursin("/.", r) || continue
        push!(sources, (relpath(joinpath(r, f), root), read(joinpath(r, f), String)))
    end
    for (r, _, fs) in walkdir(joinpath(root, "src")), f in fs
        endswith(f, ".jl") || continue
        text = read(joinpath(r, f), String)
        # A docstring of this package is a non-raw string: its LaTeX backslashes
        # are doubled in the source.
        docstrings = [replace(m.captures[1], "\\\\" => "\\") for m in eachmatch(r"\"\"\"(.*?)\"\"\""s, text)]
        push!(sources, (relpath(joinpath(r, f), root), join(docstrings, "\n")))
    end
    for (path, text) in sources, f in formulas(text)
        bad = [why for (rx, why) in refused if occursin(rx, f)]
        isempty(bad) || @error "a formula of $path breaks the notation" f bad
        @test isempty(bad)
    end
end
