using LinearAlgebra
using JSON
using DelimitedFiles
using Dierckx

struct AppDataStruct
    name::String
    nb::Int
    ngc::Int
    nh::Int
    nc::Int
    NTSDA::Int
    SJDT::Array{Any}
    SMDT::Array
    STSDAT
    q0::Array
    qd0::Array
end

function fit_xycurve(csvfile="",degree=3)
    # Define the function to get the directory of the file
    file_dir = @__DIR__ # Implement this function as needed

    # Load data from CSV
    full_path = joinpath(file_dir, csvfile)
    arr = readdlm(full_path, ',', skipstart=1)
    x = arr[:, 1]/1000.0
    y = arr[:, 2]/1000.0
    #z = arr[:, 3]
    sort_indices = sortperm(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    # Fit the spline
    spline = Spline1D(x_sorted, y_sorted, k=degree)
    return spline
end

function fit_xycurve_PUSH(csvfile="",degree=3)
    # Define the function to get the directory of the file
    file_dir = @__DIR__ # Implement this function as needed

    # Load data from CSV
    full_path = joinpath(file_dir, csvfile)
    arr = readdlm(full_path, ',', skipstart=1)
    x = arr[:, 1]
    y = arr[:, 2]
    #z = arr[:, 3]
    sort_indices = sortperm(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    # Fit the spline
    spline = Spline1D(x_sorted, y_sorted, k=degree)
    return spline
end

include("AD_contact_door310.jl")
include("AD_contact_door320.jl")
include("AD_contact_door330.jl")
include("AD_contact_door340.jl")
include("AD_220.jl")

function appdata(app)
    if app >= 310 && app < 320
        return AD310(app)
    end
    if app >= 320 && app < 329
        return AD320(app)
    end
    if app >= 330 && app < 339
        return AD330(app)
    end
    if app >= 340 && app < 349
        return AD340(app)
    end
    if app >= 220 && app < 229
        return AD220(app)
    end
end