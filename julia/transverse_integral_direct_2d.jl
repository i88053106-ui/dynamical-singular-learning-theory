using CSV
using DataFrames
using HCubature
using Printf
using QuadGK
using SpecialFunctions

const TARGET = sqrt(pi) / 3

# Success threshold for the relative difference between the independent
# two-dimensional integral and the one-dimensional Bessel representation.
const PASS_TOL = 1e-7

"""
    direct_2d(T; lo, hi, rtol, atol)

Independent evaluation of

    I(T) = ∫∫_{R^2} exp(-(y^6 + y^2 z^2 + z^6)/T) dy dz

with **no reference** to the one-dimensional Bessel representation. The
fourfold symmetry restricts the integral to the first quadrant, and the
logarithmic map

    y = exp(u),  z = exp(v),   dy dz = exp(u + v) du dv

sends (0, ∞)^2 to a finite box [lo, hi]^2. The integrand is sharply
peaked, so the box is split into panels (finer near u, v ≈ [-3, 0]) to
resolve the peak from the start rather than relying solely on adaptive
subdivision from the box centre.

The truncated tails satisfy exp(u + v) ≤ exp(2 lo), which is far below
`atol` for the default lo, so their omission is negligible.
"""
function direct_2d(
    T::Float64;
    lo::Float64=-30.0,
    hi::Float64=2.0,
    rtol::Float64=1e-9,
    atol::Float64=1e-14,
)
    T > 0 || throw(ArgumentError("T must be positive"))

    function integrand(p)
        u = p[1]
        v = p[2]
        # exp(6u) underflows smoothly to 0 for very negative u (never Inf),
        # so no 0 * Inf can occur here.
        K = exp(6 * u) + exp(2 * u + 2 * v) + exp(6 * v)

        return 4.0 * exp(u + v - K / T)
    end

    breaks = sort(unique(clamp.(
        Float64[lo, -8.0, -4.0, -2.0, -1.0, 0.0, hi],
        lo,
        hi,
    )))

    total = 0.0
    total_error = 0.0

    for i in 1:(length(breaks) - 1)
        for j in 1:(length(breaks) - 1)
            a = (breaks[i], breaks[j])
            b = (breaks[i + 1], breaks[j + 1])

            value, error = hcubature(
                integrand,
                a,
                b;
                rtol=rtol,
                atol=atol,
                maxevals=10_000_000,
            )

            total += value
            total_error += error
        end
    end

    return total, total_error
end

"""
    bessel_1d(T)

Independent one-dimensional value used only for comparison. This is a
separate computation of the Bessel representation

    I(T) = (4/3) sqrt(T) ∫ exp(-x^2) K₀(2 sqrt(T) x^3) dx;

its value is never fed into `direct_2d`.
"""
function bessel_1d(T::Float64; rtol::Float64=1e-11, atol::Float64=1e-14)
    function integrand(s)
        x = exp(s)
        return exp(s - x^2) * besselk(0, 2 * sqrt(T) * x^3)
    end

    value, _ = quadgk(integrand, -40.0, 6.0; rtol=rtol, atol=atol)

    return (4 / 3) * sqrt(T) * value
end

function main()
    temperatures = [1e-2, 1e-3, 1e-4]

    rows = DataFrame(
        T=Float64[],
        integral_direct=Float64[],
        estimated_error=Float64[],
        integral_bessel=Float64[],
        relative_difference=Float64[],
        pass=Bool[],
    )

    all_pass = true

    for T in temperatures
        direct, error = direct_2d(T)
        bessel = bessel_1d(T)

        relative_difference = abs(direct - bessel) / abs(bessel)
        passed = relative_difference < PASS_TOL
        all_pass &= passed

        push!(
            rows,
            (T, direct, error, bessel, relative_difference, passed),
        )

        @printf(
            "T=%8.1e  direct=%18.12e  bessel=%18.12e  rel=%10.3e  %s\n",
            T,
            direct,
            bessel,
            relative_difference,
            passed ? "PASS" : "FAIL",
        )
    end

    mkpath("data")
    CSV.write("data/transverse_integral_direct_2d.csv", rows)

    println()
    println("Success threshold PASS_TOL = ", PASS_TOL)
    println(all_pass ? "ALL PASS" : "SOME FAILED")
    println("Saved data/transverse_integral_direct_2d.csv")

    return all_pass
end

main()
