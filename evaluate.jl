A = [i/10 for i in 5:5:100]
include("evaluate_base.jl")
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
    simu_label = "init"
end

if exp_label== "general" # self-citing
    run(`julia evaluate_1.jl wt_15mM_salt $simu_label`)
    run(`julia evaluate_1.jl wt_150mM_salt $simu_label`)
else
    initialize(exp_label,simu_label)
    length_of_paras = length(index_p.k_on)
    @showprogress 1 "analyzing $exp_label..." for i in 1:length_of_paras
        local k_on,k_off,v_open,v_close = index_p[i,:]
        evaluate(df,k_on,k_off,v_open,v_close,folds)
        # @printf("analyzing progress: %.2f\r",i/length_of_paras)
        
    end
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
    println(@sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β))
    # k_on,k_off,v_open,v_close,α,β= 1e-5,1e-4,1e-2,1e-3,1,2
    ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β)
    yaxis!(:flip)
    ylims!(-2.0,2.0)
    xlims!(1500,2400)
    savefig("./figs/rsa_state_transition_$(exp_label)_$simu_label.svg")
    microstates_plot(k_on,k_off,v_open,v_close,folds)

## work completed
CSV.write("./data_simu/diff_$(exp_label)_$simu_label.csv",df)
end