print("initializing...\r")
using DataFrames
using ProgressMeter
using CSV
para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
# push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])

CSV.write("figs/sources/para.csv",para_df)
N = 100
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"

exp_folds = [0,1,4,10,25,50]
simu_folds = [0,1,4,10,25,50]

simu_label = rand([i for i in 1:1000])+4000


# select L of interest:
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
print("initialized.....\n")
for exp_label in exp_labels
    L = 5000
    print("processing $exp_label...\n")
    requested_df=para_df[[para_df.L[i] == L && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
    k_on,k_off,k_open,k_close,α,β,L,exp_label,loss = requested_df[1,:]
    for fold in simu_folds
        cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label) $gaps_type`
        run(cmd)
    end
    print("simulation for $exp_label completed.\n")
    include("evaluate_base.jl")
    initialize(exp_label,simu_label)
    print("\n$exp_label data loaded...\n")
    ensemble_plot(k_on,k_off,k_open,k_close,simu_folds,exp_folds,α,β,exp_label,simu_label)
    yaxis!(:flip)
    ylims!(-2.0,2.0)
    xlims!(1500,2400)
    savefig("./figs/predict_$(exp_label)_$simu_label.svg")
    microstates_plot(k_on,k_off,k_open,k_close,simu_folds,exp_label,simu_label)
    print("$exp_label plot completed.\n")
    run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`)
end
println("tasks completed.")