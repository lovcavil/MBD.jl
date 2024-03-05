using MacroTools
using Parameters
export myFunction, myConstant  # 导出任何需要的函数或常量

# Define a struct for simulation parameters
struct SimulationParams
    m1::Float64
    g::Float64
    dt::Float64
    t::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    du₀::Vector{Float64}
    u₀:: Vector{Float64}
    differential_vars::Union{Vector{Bool}, Nothing}
    equation::Union{Function, Nothing}
    SMDT
    STSDAT
    SJDT
    par
end
"""
    length = 10
    ranges = [(2, 4), (6, 8)]
    result = create_bool_vector(length, ranges)
"""
function create_bool_vector(nb::Int, nc::Int)
    differential_vars = vcat(repeat([true], 7*nb), repeat([false], nc), repeat([true], 7*nb))
    return differential_vars
end

