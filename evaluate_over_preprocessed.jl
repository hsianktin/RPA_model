A = [i/10 for i in 5:5:100]
include("preprocessed_evaluate_base.jl")
using ProgressMeter

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
    exp_label = "wt_15mM_salt"
    simu_label = "update_1"
end

if exp_label== "general" # self-citing
    run(`julia evaluate_1.jl wt_15mM_salt $simu_label`)
    run(`julia evaluate_1.jl wt_150mM_salt $simu_label`)
else
    initialize(exp_label,simu_label)
    length_of_paras = length(index_p.k_on)
    norm_debug = false
    @showprogress 1 "analyzing $exp_label..." for i in 1:length_of_paras
        local k_on,k_off,v_open,v_close = index_p[i,:]
        evaluate(df,k_on,k_off,v_open,v_close,folds)
        # @printf("analyzing progress: %.2f\r",i/length_of_paras)
    end
    
    ## debugging..
    # i = 1
    # α = A[i]
    # β = 2α 
    # k_on,k_off,v_open,v_close = index_p[i,:]
    # fold = folds[i]
    # paras = [k_on,k_off,v_open,v_close,fold,α,β]
    # paras = [convert(Float64,x) for x in paras]
    # t_X,μ_X,σ_X = access_trace_statistics(paras)
    # plot(t_X,μ_X,label="simu_$fold")
    # i=2
    # fold = folds[i]
    # paras = [k_on,k_off,v_open,v_close,fold,α,β]
    # paras = [convert(Float64,x) for x in paras]
    # t_X,μ_X,σ_X = access_trace_statistics(paras)
    # plot!(t_X,μ_X,label="simu_$fold")

    # test_df = CSV.read("./data/simu_$(exp_label)_$(simu_label)_$(k_on)_$(k_off)_$(v_open)_$(v_close).csv",DataFrame)
    # simu_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_1 = Array{Float64,1}[],σ_1 = Array{Float64,1}[], μ_2 = Array{Float64,1}[],σ_2 = Array{Float64,1}[])
    # @show test_df
    # for i in 1:length(test_df.folds)
    #     fold′ = test_df.folds[i]
    #     T = [parse(Float64,x) for x in split(test_df.T[i][2:end-1],",")]
    #     μ_1 = [parse(Float64,x) for x in split(test_df.μ_1[i][2:end-1],",")] 
    #     σ_1 = [parse(Float64,x) for x in split(test_df.σ_1[i][2:end-1],",")]
    #     μ_2 = [parse(Float64,x) for x in split(test_df.μ_2[i][2:end-1],",")]
    #     σ_2 = [parse(Float64,x) for x in split(test_df.σ_2[i][2:end-1],",")]
    #     push!(simu_df, [fold′, T, μ_1, σ_1, μ_2, σ_2])   
    # end
    # @show simu_df
    # j = findfirst(x-> x== fold, simu_df.folds)
    # T = simu_df.T[j]
    # μ_1 = simu_df.μ_1[j]
    # σ_1 = simu_df.σ_1[j]
    # μ_2 = simu_df.μ_2[j]
    # σ_2 = simu_df.σ_2[j]
    # μ_X = α.*μ_1 + β.*μ_2
    # σ_X = sqrt.(α.*(σ_1.^2) + β.*(σ_2.^2))
    # end of debugging
    # diff contains all the information needed
    if length(unique(df.k_on)) > 1
        analyze(df,"k_on")
    end
    if length(unique(df.k_off)) > 1
        analyze(df,"k_off")
    end
    if length(unique(df.v_open)) > 1
        analyze(df,"v_open")
    end
    if length(unique(df.v_close)) > 1
        analyze(df,"v_close")
    end
    analyze(df,"α")
    analyze(df,"β")
    # CSV.write("./data_simu/analyze_$(exp_label)_$(simu_label).csv")
    minval,index_min = findmin(df.diff)
    println("minval = $minval")
    k_on,k_off,v_open,v_close,α,β,d=df[index_min,:]
    norm_debug = true
    diff(k_on,k_off,v_open,v_close,folds,α,β)
    println(@sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β))
    # k_on,k_off,v_open,v_close,α,β= 1e-5,1e-4,1e-2,1e-3,1,2
    ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β)
    yaxis!(:flip)
    ylims!(0.5,2.0)
    xlims!(1500,2400)
    savefig("./figs/rsa_state_transition_$(exp_label)_$simu_label.svg")
    microstates_plot(k_on,k_off,v_open,v_close,folds)
end