using CSV
using DataFrames
using Printf
using QuadGK
using SpecialFunctions

const TARGET = sqrt(pi) / 3
const EULER_GAMMA = Base.MathConstants.eulergamma

"""
    transverse_integral(T)

Compute

I(T) = ∫∫ exp(-(y^6 + y^2 z^2 + z^6)/T) dy dz

using the exact one-dimensional representation

I(T) = (4/3) sqrt(T) ∫ exp(-x^2) K₀(2 sqrt(T) x^3) dx.
"""
function transverse_integral(T::Float64)
    T > 0 || throw(ArgumentError("T must be positive"))

    # x = exp(s) removes the logarithmic endpoint singularity at x = 0.
    function integrand(s)
        x = exp(s)
        argument = 2 * sqrt(T) * x^3

        return exp(s - x^2) * besselk(0, argument)
    end

    # The omitted tails are negligible at Float64 precision.
    value, error = quadgk(
        integrand,
        -40.0,
        6.0;
        rtol=1e-11,
        atol=1e-14,
    )

    scale = (4 / 3) * sqrt(T)

    return scale * value, scale * error
end

function main()
    exponents = range(2.0, 12.0; length=21)
    temperatures = 10.0 .^ (-collect(exponents))

    rows = DataFrame(
        T=Float64[],
        integral=Float64[],
        estimated_error=Float64[],
        raw_ratio=Float64[],
        corrected_ratio=Float64[],
        target=Float64[],
    )

    correction_constant = EULER_GAMMA + 6 * log(2)

    for T in temperatures
        value, error = transverse_integral(T)

        raw_denominator =
            sqrt(T) * log(1 / T)

        corrected_denominator =
            sqrt(T) * (log(1 / T) + correction_constant)

        raw_ratio = value / raw_denominator
        corrected_ratio = value / corrected_denominator

        push!(
            rows,
            (
                T,
                value,
                error,
                raw_ratio,
                corrected_ratio,
                TARGET,
            ),
        )

        @printf(
            "T=%10.3e  I(T)=%14.7e  raw=%10.7f  corrected=%10.7f\n",
            T,
            value,
            raw_ratio,
            corrected_ratio,
        )
    end

    mkpath("data")
    CSV.write("data/transverse_integral_julia.csv", rows)

    println()
    println("Target sqrt(pi)/3 = ", TARGET)
    println("Saved data/transverse_integral_julia.csv")
end

main()
