# initialize necessary environment
println("initializing modules...")
push!(LOAD_PATH,pwd())
using TonksGaswithReactionGaps
try
    using Statistics
    using Random
    using CSV
    using Plots
    using DataFrames
    using Printf
    using ProgressMeter
    using Distributed
    using DelimitedFiles
catch
    println("dependency requirements are not met. Try installing packages...")
    using Pkg
    Pkg.add("Statistics")
    Pkg.add("CSV")
    Pkg.add("Plots")
    Pkg.add("DataFrames")
    Pkg.add("ProgressMeter")
    Pkg.add("DelimitedFiles")
    using Statistics
    using Random
    using CSV
    using Plots
    using DataFrames
    using Printf
    using ProgressMeter
    using Distributed
    using DelimitedFiles
end
println("module initialized")
println("initializing folders...")
try 
    mkdir("data_simu")
catch
    println("data_simu/ already exists")
end

try
    mkdir("data_exp")
catch
    println("data_exp/ already exists")
end

try
    mkdir("figs")
    mkdir("figs/landscape")
    mkdir("figs/sources")
catch
    println("figs/ already exists")
end

println("initialization completed")
