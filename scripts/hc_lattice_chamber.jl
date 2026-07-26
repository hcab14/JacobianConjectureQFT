# Q2a of docs/CLASSICAL_MAP_INVARIANTS.md §5.2: real chamber function of the
# 2-site lattice deformation of the Alpöge–Mathew map,
#
#   F_eps(phi)_0 = F(phi_0) + eps (phi_1 - phi_0)
#   F_eps(phi)_1 = F(phi_1) + eps (phi_0 - phi_1)
#
# (convention of scripts/classical_map_invariants_probe.py: F_M^K = F_M + K.phi
# with K = eps [[-I, I], [I, -I]]).  Driven by scripts/lattice_chamber.py,
# which generates the input file (exact eps=0 product start points from
# jcqft.fibers.exact_fiber, targets, eps ladder) and parses the output.
#
# Usage:  julia scripts/hc_lattice_chamber.jl /tmp/lattice_chamber/input.jl
#
# METHOD.  Three mutually cross-checking computations per (target J, eps):
#  (m) MASTER: treat (eps, J) in C^7 as parameters; polyhedral solve at a
#      generic complex point + monodromy stabilization gives the master
#      solution set (D paths); every fixed (eps, J) is then reached by
#      parameter homotopy of ALL D paths — this sees solutions that come in
#      from infinity relative to eps=0 (the interesting phenomenon).
#  (a) eps-HOMOTOPY: the 9 exact product solutions at eps=0 are tracked to
#      (eps, J) through a complex detour; every finite endpoint must
#      reappear in (m) — and the difference (m) minus (a) counts solutions
#      that entered from infinity.
#  (b) FRESH polyhedral solve at the fixed rational (eps, J); any solution
#      NOT already in (m) would prove (m) incomplete -> hard failure.
# Reality/distinctness is CERTIFIED (HC.certify, interval arithmetic) against
# the exact rational system.  Completeness of (m) is NOT certified (standard
# homotopy-continuation caveat); it is cross-checked by (b) and by transport
# of the master set to an independent generic point.
#
# Output: machine-readable MASTER / STARTS / RESULT / BISECT lines (tab-
# separated), plus [ok]/[FAIL] lines.  Exits nonzero on any inconsistency.

using HomotopyContinuation
using LinearAlgebra
using Random

if length(ARGS) < 1
    println("usage: julia scripts/hc_lattice_chamber.jl <input.jl>")
    exit(2)
end
include(abspath(ARGS[1]))   # defines eps_list, targets, bisect

const T0 = time()
NFAIL = 0

function check(label, cond)
    global NFAIL
    if cond
        println("  [ok] $label   ($(round(time() - T0, digits=1)) s)")
    else
        NFAIL += 1
        println("  [FAIL] $label   ($(round(time() - T0, digits=1)) s)")
    end
    cond
end

# --- the lattice map --------------------------------------------------------
@var x0 y0 z0 x1 y1 z1
@var e a0 b0 c0 a1 b1 c1

AM(x, y, z) = [
    (1 + x*y)^3*z + y^2*(1 + x*y)*(4 + 3*x*y),
    y + 3*x*(1 + x*y)^2*z + 3*x*y^2*(4 + 3*x*y),
    2*x - 3*x^2*y - x^3*z,
]

vars = [x0, y0, z0, x1, y1, z1]
phi0 = [x0, y0, z0]
phi1 = [x1, y1, z1]
pars = [e, a0, b0, c0, a1, b1, c1]

eqs = vcat(
    AM(x0, y0, z0) .+ e .* (phi1 .- phi0) .- [a0, b0, c0],
    AM(x1, y1, z1) .+ e .* (phi0 .- phi1) .- [a1, b1, c1],
)
P = System(eqs; variables = vars, parameters = pars)

matches(s, B; tol = 1e-6) = any(b -> norm(s - b) < tol * max(1.0, norm(s)), B)

# --- 0. verify the exact eps=0 product start points -------------------------
println("=== 0. eps=0 start points (exact product fibers from Python) ===")
for tg in targets
    par0 = ComplexF64[0; tg.J]
    maxres = maximum(s -> maximum(abs.(P(s, par0))), tg.starts)
    check("$(tg.name): 9 product starts, max eps=0 residual = " *
          "$(round(maxres, sigdigits=2)) < 1e-8",
          length(tg.starts) == 9 && maxres < 1e-8)
    nreal0 = count(s -> maximum(abs.(imag.(s))) < 1e-10, tg.starts)
    check("$(tg.name): eps=0 real product count = $nreal0 = expected $(tg.N0)",
          nreal0 == tg.N0)
    println("STARTS\t$(tg.name)\t9\t$nreal0")
end

# --- 1. master solution set over generic complex (eps, J) -------------------
println("=== 1. master set: monodromy over the joint (eps, J) parameter space ===")
Random.seed!(20260726)
p0 = randn(ComplexF64, 7)
r0 = solve(P; target_parameters = p0, start_system = :polyhedral,
           show_progress = false)
mr = monodromy_solve(P, solutions(r0), p0; show_progress = false,
                     max_loops_no_progress = 20)
S0 = solutions(mr)
D = length(S0)
println("MASTER\tD=$D\tpolyhedral_seed=$(nsolutions(r0))")
check("master degree D = $D >= polyhedral seed $(nsolutions(r0)) " *
      "(monodromy found the extra paths)", D >= nsolutions(r0))
p1 = randn(ComplexF64, 7)
rt = solve(P, S0; start_parameters = p0, target_parameters = p1,
           show_progress = false)
check("master set transports to an independent generic point without loss " *
      "($(nsolutions(rt)) of $D)", nsolutions(rt) == D)

# --- helper: certified fiber at exact rational (eps, J) ---------------------
# Three independent transport routes of the master set (direct, and via two
# random complex midpoints) are unioned before certification, so a lost
# path on one route is recovered by the others; certify dedupes.  The RNG
# is reseeded deterministically per point so counts are reproducible
# across runs regardless of how many points were evaluated before.
function fiber_at(epsr::Rational{Int}, Jr::Vector{Rational{Int}})
    Random.seed!(hash((epsr, Jr)))
    pt = ComplexF64[epsr; Jr]
    rA = solve(P, S0; start_parameters = p0, target_parameters = pt,
               show_progress = false)
    sols = solutions(rA; only_nonsingular = false)
    for _ in 1:2
        pmid = randn(ComplexF64, 7)
        rM = solve(P, S0; start_parameters = p0, target_parameters = pmid,
                   show_progress = false)
        rB = solve(P, solutions(rM); start_parameters = pmid,
                   target_parameters = pt, show_progress = false)
        append!(sols, solutions(rB; only_nonsingular = false))
    end
    cert = certify(P, sols, [epsr; Jr])       # exact rational certification
    (sols = sols,
     n_complex = ndistinct_certified(cert),
     n_real = ndistinct_real_certified(cert))
end

# --- 2. results per (target, eps) -------------------------------------------
println("=== 2. chamber counts per (target, eps): master + eps-homotopy + fresh ===")
detour = 0.62 + 0.79im    # complex detour for the eps-homotopy legs
RES = Dict{Tuple{String,Rational{Int}},Tuple{Int,Int}}()
for tg in targets
    for epsr in eps_list
        fib = fiber_at(epsr, tg.J)
        RES[(tg.name, epsr)] = (fib.n_complex, fib.n_real)

        # (a) eps-homotopy from the 9 exact eps=0 product solutions
        mid = ComplexF64[detour * epsr; tg.J]
        fin = ComplexF64[epsr; tg.J]
        ra = solve(P, tg.starts; start_parameters = ComplexF64[0; tg.J],
                   target_parameters = mid, show_progress = false)
        rb = solve(P, solutions(ra); start_parameters = mid,
                   target_parameters = fin, show_progress = false)
        eps0_ends = solutions(rb; only_nonsingular = false)
        n_eps0 = length(eps0_ends)
        n_eps0_matched = count(s -> matches(s, fib.sols), eps0_ends)
        n_eps0_real = count(
            s -> maximum(abs.(imag.(s))) < 1e-8 && matches(s, fib.sols),
            eps0_ends)

        # (b) fresh polyhedral solve at the fixed rational (eps, J)
        Pfix = System(subs(eqs, pars => [epsr; Rational{Int}.(tg.J)]);
                      variables = vars)
        rf = solve(Pfix; start_system = :polyhedral, show_progress = false)
        fresh = [s for s in solutions(rf; only_nonsingular = false)
                 if maximum(abs.(Pfix(s))) < 1e-6 && norm(s) < 1e8]
        n_fresh_new = count(s -> !matches(s, fib.sols), fresh)

        check("$(tg.name) eps=$epsr: all $n_eps0 finite eps-homotopy " *
              "endpoints reappear in the master fiber",
              n_eps0_matched == n_eps0)
        check("$(tg.name) eps=$epsr: fresh polyhedral solve " *
              "($(length(fresh)) sols) finds nothing outside the master " *
              "fiber", n_fresh_new == 0)
        check("$(tg.name) eps=$epsr: conjugation parity " *
              "(n_complex - n_real = $(fib.n_complex - fib.n_real) even)",
              iseven(fib.n_complex - fib.n_real))
        println("RESULT\t$(tg.name)\t$epsr\t$D\t$(fib.n_complex)\t" *
                "$(fib.n_real)\t$n_eps0\t$n_eps0_real\t$(length(fresh))")
    end
end

# --- 3. bisection of a real-count jump in J at fixed eps ---------------------
println("=== 3. bisection: wall location on the segment " *
        "$(bisect.from) -> $(bisect.to) at eps = $(bisect.eps) ===")
Jlo = Rational{Int}.(targets[findfirst(t -> t.name == bisect.from, targets)].J)
Jhi = Rational{Int}.(targets[findfirst(t -> t.name == bisect.to, targets)].J)
Jseg(t) = (1 - t) .* Jlo .+ t .* Jhi
Nlo = fiber_at(bisect.eps, Jseg(0//1)).n_real
Nhi = fiber_at(bisect.eps, Jseg(1//1)).n_real
check("bisection endpoints have different certified real counts " *
      "($Nlo vs $Nhi)", Nlo != Nhi)
lo, hi = 0//1, 1//1
while hi - lo > bisect.tol
    tm = (lo + hi) // 2
    Nm = fiber_at(bisect.eps, Jseg(tm)).n_real
    println("  bisect t=$tm: N_real = $Nm")
    if Nm == Nlo
        global lo = tm
    else
        global hi = tm
    end
end
println("BISECT\t$(bisect.eps)\t$(bisect.from)\t$(bisect.to)\t$lo\t$hi\t$Nlo\t$Nhi")
println("  wall crossing bracketed in t in ($lo, $hi), " *
        "N_real: $Nlo -> $Nhi")

# --- 4. summary --------------------------------------------------------------
println("=== 4. summary table (certified real counts; complex in parens) ===")
hdr = "  eps        " * join([lpad(t.name, 12) for t in targets])
println(hdr)
println("  0 (exact)  " *
        join([lpad("$(t.N0)  (9 C)", 12) for t in targets]))
for epsr in eps_list
    row = "  " * rpad("$epsr", 11)
    for t in targets
        nc, nr = RES[(t.name, epsr)]
        row *= lpad("$nr ($nc C)", 12)
    end
    println(row)
end

if NFAIL > 0
    println("hc_lattice_chamber: $NFAIL INTERNAL INCONSISTENCIES")
    exit(1)
end
println("hc_lattice_chamber: all internal cross-checks passed " *
        "($(round(time() - T0, digits=1)) s)")
