# the original evaluation script is considerably slow,
# this script is to speed up the evaluation process by
# using the split-apply-combine strategy and multi-threads parallelization

using Statistics
using Pipe
using CSV
using Plots
using DataFrames
using LinearAlgebra
using DelimitedFiles
using Printf
∑(x) = sum(x)
include("evaluate_base.jl")
#### for test purpose
simupath = "$(pwd())/data_simu"
# data_path = raw""
T = 2400
exp_label = "wt_15mM_salt"
simu_label = "update_1"
####
if length(ARGS) >= 3
    exp_label = ARGS[1]
    simu_label = ARGS[2]
    T = parse(Float64, ARGS[3])
else
    println("Usage: julia fit_base.jl <exp_label> <simu_label> <T>")
end
println("exp_label = $exp_label")
println("simu_label = $simu_label")
println("T = $T")
println("Threads: $(Threads.nthreads())")
function get_dbpath(exp_label,simu_label)
    path = "$simupath/rsa_plot_$(exp_label)_$(simu_label).csv"
    return path
end

# need to specify total time T, assume that the data is stored per minute
function load_simu_df(dbpath,T)
    # rewrite the load function using  DataFrames
    para_headers = [:k_on,:k_off,:v_open,:v_close,:fold,:L,:T1,:T2,:N,:mode]
    headers = vcat(para_headers, [Symbol(i) for i ∈ 0:T])
    df = CSV.read(dbpath, header=headers, DataFrame)
    ids = Array{Int}(undef,nrow(df))
    id = 0
    init_fold = df.fold[1]
    for i in 1:nrow(df)
        if i % 2 == 1
            if df.fold[i] == init_fold
                id += 1
            else
                id = 1
                init_fold = df.fold[i]
            end
        end
        ids[i] = id
    end
    df.id = ids
    mode(x) = x == 0 ? "fraction_30nt" : "fraction_20nt"
    df.modes .= mode.(df.mode)
    para_headers = [:id,:k_on,:k_off,:v_open,:v_close,:fold,:L,:T1,:T2,:N,:modes]
    compact_para_headers = [:id,:k_on,:k_off,:v_open,:v_close,:fold,:L,:T1,:T2,:N, :time]
    long_df = @pipe df |>
                stack(
                    _, 
                    [Symbol(i) for i ∈ 0:T], 
                    para_headers, 
                    variable_name=:time, 
                    value_name=:fraction,
                ) |>
                transform(
                    _, 
                    :time => ByRow(x -> parse(Float64,x)) => :time
                    ) |>
                unstack(
                     _,
                     compact_para_headers,
                     :modes,
                     :fraction,
                     allowduplicates=false
                    )
    compact_para_headers = [:k_on,:k_off,:v_open,:v_close,:fold,:L,:T1,:T2,:N, :time]
    microstates_stats = @pipe long_df |> 
        groupby(_, compact_para_headers) |>
        combine(_, 
               :fraction_30nt => mean =>:f̄₃₀,
               :fraction_30nt => std => :Δf₃₀,
               :fraction_20nt => mean =>:f̄₂₀,
               :fraction_20nt => std => :Δf₂₀,
              ) |>
        sort(
            _,
            compact_para_headers
        )
    return microstates_stats
    # the output is a DataFrame with the following columns:
    # k_on, k_off, v_open, v_close, fold, L, T1, T2, N, time, f̄₃₀, Δf₃₀, f̄₂₀, Δf₂₀
end

function load_exp_df(folds, exp_dict,T̄)
    exp_df = DataFrame(
        fold = Float64[], 
        ℓ̄ = Float64[], 
        Δℓ̄ = Float64[], 
        Ī = Float64[], 
        time = Float64[]
        )
    for fold ∈ folds
        conc = convert(Float64, fold)
        T, L , ΔL = access_trace_statistics(exp_dict["$conc"][1])
        I, T_I = access_intensity_trace(exp_label, fold)
        for i ∈ eachindex(T)
            if T[i] ≤ T̄
                push!(exp_df, [fold, L[i], ΔL[i], I[i] , T[i]])
            end
        end
    end
    return exp_df
end



function initialize(exp_label,dbpath,T)
    global exp_dict = Dict{String,Array{String,1}}()
    if exp_label == "wt_15mM_salt"
        global exp_data_base=[("wt",0),("wt",1),("wt",4),("wt",10),("wt",25),("wt",50)]
        global folds = [0,1,4,10,25]
    elseif exp_label == "wt_150mM_salt"
        global exp_data_base=[("lwt",0),("lwt",1),("lwt",4),("lwt",10),("lwt",25),("lwt",50)]
        global folds = [0,1,4,10,25]
    end
    function salt_concentration(type) # historically, the wt label refers to wt at 15mM salt; lwt for 150mM salt.
        if type == "wt"
            return "15mM salt"
        elseif type == "lwt"
            return "150mM salt"
        end
    end
    function exp_dict_inject(type::String,concentration::Float64)
        datapath = "./data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$concentration"]=[datapath,fig_label]
    end
    function exp_dict_inject(type::String,concentration::Int)
        datapath = "./data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$(convert(Float64,concentration))"]=[datapath,fig_label]
    end
    function exp_dict_inject(tuple::Tuple)
        if length(tuple) == 2
            type,concentration = tuple
            exp_dict_inject(type,concentration)
        else
            println("input not match")
        end
    end
    exp_dict_inject.(exp_data_base)
    # dbpath = dbpath(exp_label,simu_label)
    simu_df = load_simu_df(dbpath,T)
    exp_df = load_exp_df(folds, exp_dict, T)
    return simu_df, exp_df
end

function fit(folds,f₃₀,f₂₀,ℓ,I,times)
    local df = DataFrame(fold=folds, f₃₀=f₃₀, f₂₀=f₂₀, ℓ=ℓ, I = I,  time=times)
    loss_df = DataFrame(α=Float64[], β=Float64[], difference=Float64[])
    A = [i for i in 0.5:0.1:4]
    for α in A
        B = [i*α/5 for i in 7:20]
        for β in B
        # β = 2*α
            difference = 0
            for fold in unique(df.fold)
                df_slice = @pipe df |>
                    filter(row -> row.fold == fold && row.time ≥ 1200, _)
                y = df_slice.ℓ
                i = df_slice.I
                f₃₀ = df_slice.f₃₀
                f₂₀ = df_slice.f₂₀
                t = df_slice.time
                ŷ = (α*f₃₀)+(β*f₂₀)+(-f₃₀-f₂₀.+1) |> (x -> x./mean([x[convert(Int,1775/5+1-1200/5)],x[convert(Int,1790/5+1-1200/5)]]))
                Î = f₃₀ + f₂₀ |> (x -> x./mean([x[convert(Int,1775/5+1-1200/5)],x[convert(Int,1790/5+1-1200/5)]]))
                εₗ = y-ŷ
                εᵢ = i-Î
                difference += norm(εₗ)^2 + norm(εₗ[t .≥ 2000])^2 + (norm(εᵢ)^2 + norm(εᵢ[t .≥ 2000])^2) * intensity_weight^2
            end
            push!(loss_df,[α,β,difference])
            # println([k_on,k_off,v_open,v_close,α,β,difference])
        end
        # only enabled in lab 
        # ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β);
        # savefig("D:\\rpa\\figs\\plot_$(k_on)_$(k_off)_$(v_open)_$(v_close)_$(α)_$(β).png")
    end
    α,β,difference = @pipe loss_df |> sort(_,:difference) |> x -> x[1,:]
    return (α,β,difference)
end

function analyze_df(microstates_stats, exp_df)
    times = exp_df.time
    simu_df = filter(row -> row.time ∈ times, microstates_stats)
    # by doing so, the simu_df are sampled at the same time points as the exp_df
    
    # require that the microstates_stats are well sorted
    ℓ = exp_df.ℓ̄
    I = exp_df.Ī
    para_headers = [:k_on,:k_off,:v_open,:v_close,:L,:T1,:T2,:N]
    combine_df = @pipe simu_df |>
        groupby(_, para_headers) |>
        combine(_,
            [:fold, :f̄₃₀, :f̄₂₀] => ((x,y,z) -> (α=fit(x,y,z,ℓ,I,times)[1],β=fit(x,y,z,ℓ,I,times)[2],loss=fit(x,y,z,ℓ,I,times)[3])) => AsTable,
        )
    result = @pipe combine_df |>
            sort(_,:loss) |>
            select(_,[:k_on,:k_off,:v_open,:v_close,:α,:β,:loss]) |>
            x -> x[1,:]
    return result
end
function visualize(result,simu_df)
    k_on, k_off, v_open, v_close, α, β, loss = result[1,:]
    combine_df = @pipe simu_df |> # find the entries that are defined by the result df
        filter(row -> row.k_on == k_on && row.k_off == k_off && row.v_open == v_open && row.v_close == v_close , _) |>
        transform(_, [:N] => (x -> (α=α, β=β, loss=loss)) => AsTable) |> 
        transform(_, 
            [:f̄₃₀, :f̄₂₀] => ((f₃₀,f₂₀) -> (α * f₃₀ + β * f₂₀ + 1 * (-f₃₀-f₂₀.+1)) ) => :ℓ̂
            )
    # combine_df is a dataframe with ℓ̂ as the absolute length
    norm_dict = @pipe groupby(combine_df, :fold) |>
        combine(_,
            [:time, :ℓ̂] => ((t,y)-> mean(y[t .< 1800 .&& t .> 1600])) => :norm
        ) |> Dict(Pair.(_.fold,_.norm))

    plot_df = @pipe combine_df |>
    transform!(
        _,
        [:fold, :ℓ̂] => ((fold,y) -> y./[norm_dict[f] for f in fold]) => :ext
    )
    CSV.write("figs/sources/plot_df_$(simu_label)_$(exp_label).csv",plot_df)
end
data_path = get_dbpath(exp_label,simu_label)
simu_df, exp_df = initialize(exp_label,data_path,T);
result = analyze_df(simu_df, exp_df) |> DataFrame
CSV.write("figs/sources/result_$(exp_label)_$(simu_label).csv",result)
visualize(result,simu_df)
