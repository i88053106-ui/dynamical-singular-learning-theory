using QuadGK
using HCubature
using CSV
using DataFrames

expected = sqrt(pi) / 3

println("Julia version: ", VERSION)
println("sqrt(pi) / 3 = ", expected)
println("Julia environment OK")
