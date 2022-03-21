# in order to accelerate the training process, preprocess data to avoid repeated computation.
include("evaluate_base.jl")
using ProgressMeter

function process_data(exp_label,folds)
    for fold in folds
        exp_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_X = Array{Float64,1}[],σ_X = Array{Float64,1}[])
        conc = convert(Float64,fold)
        time_course, μ_X, σ_X = access_trace_statistics(exp_dict["$conc"][1])
        push!(exp_df,[fold,time_course,μ_X,σ_X])
        CSV.write("./data/exp_$(exp_label)_$(fold).csv",exp_df)
    end
end

function process_data(exp_label,simu_label,folds,k_on,k_off,v_open,v_close)
    simu_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_1 = Array{Float64,1}[],σ_1 = Array{Float64,1}[], μ_2 = Array{Float64,1}[],σ_2 = Array{Float64,1}[])
    for fold in folds
        time_course, μ_1, σ_1, μ_2, σ_2 = access_microstates(k_on,k_off,v_open,v_close,fold)
        push!(simu_df,[fold,time_course,μ_1,σ_1,μ_2,σ_2])
    end
    CSV.write("./data/simu_$(exp_label)_$(simu_label)_$(k_on)_$(k_off)_$(v_open)_$(v_close).csv",simu_df)
end

if length(ARGS) == 1
    exp_label = ""
    simu_label = ARGS[1]
    # simu_label = "wt_15mM_salt"
elseif length(ARGS) == 2
    exp_label=ARGS[1]
    println("exp_label=$exp_label")
    simu_label=ARGS[2]
    println("simu_label=$simu_label")
else
    exp_label = "wt_150mM_salt"
    simu_label = "update_1"
end

if exp_label== "general" # self-citing
    run(`julia evaluate.jl wt_15mM_salt $simu_label`)
    run(`julia evaluate.jl wt_150mM_salt $simu_label`)
else
    initialize(exp_label,simu_label)  
    println("processing experiment data for $exp_label")
    process_data(exp_label,folds)
    println("experiment data processed")
    length_of_paras = length(index_p.k_on)
    @showprogress 1 "processing $exp_label $simu_label..." for i in 1:length_of_paras
        local k_on,k_off,v_open,v_close = index_p[i,:]
        process_data(exp_label,simu_label,folds,k_on,k_off,v_open,v_close)
    end
    CSV.write("./data/index_$(exp_label)_$(simu_label).csv",index_p)
## work completed
end