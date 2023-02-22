# Based on evaluate_base, provide information about losses and landscape.
# variants include evaluate_sensitivity.jl and evaluate_diffusion.jl
using DataFrames,CSV

A = [i/10 for i in 5:5:100] # values for α to be tested

include("evaluate_base.jl")
using ProgressMeter
if length(ARGS) == 1
    exp_label = ""
    simu_label = ARGS[1]
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
    run(`julia evaluate.jl wt_15mM_salt $simu_label`)
    run(`julia evaluate.jl wt_150mM_salt $simu_label`)
else
    initialize(exp_label,simu_label)
    length_of_paras = length(index_p.k_on)
    # headless plot mode for server
    ENV["GKSwstype"] = "100"

    @showprogress 1 "analyzing $exp_label..." for i in 1:length_of_paras
        local k_on,k_off,v_open,v_close = index_p[i,:]
        evaluate(df,k_on,k_off,v_open,v_close,folds)
    end
    # diff contains all the information needed
    # if length(unique(df.k_on)) > 1
        analyze(df,"k_on")
    # end
    # if length(unique(df.k_off)) > 1
        analyze(df,"k_off")
    # end
    # if length(unique(df.v_open)) > 1
        analyze(df,"v_open")
    # end
    # if length(unique(df.v_close)) > 1
        analyze(df,"v_close")
    # end
    # if length(unique(df.α)) > 1
        analyze(df,"α")
    # end
    # if length(unique(df.β)) > 1
        analyze(df,"β")
    # end
    minval,index_min = findmin(df.diff)
    println("minval = $minval")
    k_on,k_off,v_open,v_close,α,β,d=df[index_min,:]
    println(@sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β))
    ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β)
    yaxis!(:flip)
    ylims!(-2.0,2.0)
    # xlims!(1500,2400)
    savefig("./figs/rsa_state_transition_$(exp_label)_$simu_label.svg")
    # microstates_plot(k_on,k_off,v_open,v_close,folds)
    paradf = DataFrame(
        k_on = [k_on],
        k_off = [k_off],
        v_open = [v_open],
        v_close = [v_close],
        α = [α],
        β = [β],
        loss = [minval]
    )
    f = "./figs/sources/loss_$(exp_label)_$(simu_label).csv"
    CSV.write(f,paradf)
## work completed
CSV.write("./data_simu/diff_$(exp_label)_$simu_label.csv",df)
end