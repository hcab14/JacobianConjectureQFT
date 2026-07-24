# Homotopy-continuation discovery pass over the six reduced degree-2
# (2,-1,-3) mechanism queries (docs/SEARCH_213.md §5.3).
#
# Usage:  julia scripts/hc_drive_213.jl [max_paths]
#
# For each /tmp/hc_213/query_<i>.jl: build the (overdetermined) system,
# square up, measure the start-system size first, and only track if the
# path count is below the budget (default 5_000_000).  Solutions, if
# any, are printed with residuals; certification against the ORIGINAL
# overdetermined system is attempted via `certify`.
#
# Verdict semantics: "0 solutions" here is EVIDENCE of emptiness (not a
# proof -- paths can diverge); any certified solution is a concrete
# counterexample-mechanism witness for exact follow-up.

using HomotopyContinuation

max_paths = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5_000_000

for i in 0:5
    file = "/tmp/hc_213/query_$(i).jl"
    isfile(file) || continue
    include(file)   # defines vars, polys, label
    F = System(polys; variables = vars)
    n, N = length(vars), length(polys)
    println("== query $i: $label  ($N eqs, $n unknowns) ==")
    t0 = time()
    # Square up: n random linear combinations of the N polynomials.
    G = System(randn(ComplexF64, n, N) * polys; variables = vars)
    bez = prod(degrees(G))
    println("  Bezout number of squared-up system: $bez")
    if bez > max_paths
        # Try polyhedral (mixed volume can be far below Bezout).
        println("  measuring mixed volume (polyhedral start system)...")
        mv = try
            paths_to_track(G; start_system = :polyhedral)
        catch e
            println("  mixed-volume computation failed: $(typeof(e))")
            -1
        end
        println("  polyhedral paths: $mv")
        if mv < 0 || mv > max_paths
            println("  SKIP: exceeds budget $max_paths  [$(round(time()-t0)) s]")
            println()
            continue
        end
    end
    res = solve(G; start_system = :polyhedral, show_progress = false)
    sols = solutions(res; only_nonsingular = false)
    println("  tracked: $(nresults(res)); raw solutions: $(length(sols))")
    # Filter to solutions of the FULL overdetermined system.
    good = [s for s in sols if maximum(abs.(F(s))) < 1e-8]
    println("  satisfying all $N equations (res < 1e-8): $(length(good))")
    if !isempty(good)
        cert = certify(F, good)
        println("  certified: $(ndistinct_certified(cert))")
        for s in good[1:min(end, 5)]
            println("    ", s)
        end
    end
    println("  [$(round(time()-t0)) s]")
    println()
end
println("hc_drive_213 done.")
