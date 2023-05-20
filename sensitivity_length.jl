<<<<<<< HEAD
# determine the local sensitivity with respect to lengths of DNA.
# must be run after figs/para_fitted.csv is created
print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(12) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
print("multi-threading initialized.\n")

para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
# push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])

# CSV.write("figs/sources/ini_para.csv",para_df)
N = 10 # copy per command, total 100 copies per parameter combination
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"

exp_folds = [0,1,4,10,25,50]
simu_folds = [0,1,4,10,25,50]

simupath="./data_simu"

# select L of interest:
Lengths = [500,5000,50000]
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
print("initialized.....\n")
for L in Lengths
    simu_label = "L_$(L)_$(rand([i for i in 1:1000]))"
    for exp_label in exp_labels
        cmds = Array{Cmd,1}()
        print("processing $exp_label...\n")
        requested_df=para_df[[para_df.L[i] == 5000 && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
        k_on,k_off,k_open,k_close,α,β,_,exp_label,loss = requested_df[1,:]
        count=0
        for fold in simu_folds
            for i in 1:10
                count+=1
                cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$(count) $gaps_type`
                push!(cmds,cmd)
            end
        end
        @showprogress "L=$(L) $exp_label" pmap(run,cmds)
        for i in 1:count
            # try
            open("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv","r") do input
                open("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv","a") do output
                    n = countlines(input)
                    record = 0
                    seekstart(input)
                    for k in 1:n
                        line = readline(input)
                        println(output,line)
                        record = record + 1
                    end
                    # println(record)
                end
            end
        #     catch
        #     end
        end
        for i in 1:count
            try
            rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
            catch
            end
        end
        for i in 1:count
            try
                open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv","r") do input
                    open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label).csv","a") do output
                        n = countlines(input)
                        record = 0
                        seekstart(input)
                        for k in 1:n
                            line = readline(input)
                            println(output,line)
                            record = record + 1
                        end
                        # println(record)
                    end
                end
            catch
            end
        end
        for i in 1:count
            try
                rm("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv")
            catch
            end
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
        # run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`) analysis of gaps is postponed for sake of time
    end
end
println("tasks completed.")

=======
# determine the local sensitivity with respect to lengths of DNA.
# must be run after figs/para_fitted.csv is created
print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(12) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
print("multi-threading initialized.\n")

para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
# push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])

# CSV.write("figs/sources/para.csv",para_df)
N = 10 # copy per command, total 100 copies per parameter combination
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"

exp_folds = [0,1,4,10,25,50]
simu_folds = [0,1,4,10,25,50]

simupath="./data_simu"

# select L of interest:
Lengths = [500,5000,50000]
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
print("initialized.....\n")
for L in Lengths
    simu_label = "L_$(L)_$(rand([i for i in 1:1000]))"
    for exp_label in exp_labels
        cmds = Array{Cmd,1}()
        print("processing $exp_label...\n")
        requested_df=para_df[[para_df.L[i] == 5000 && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
        k_on,k_off,k_open,k_close,α,β,_,exp_label,loss = requested_df[1,:]
        count=0
        for fold in simu_folds
            for i in 1:10
                count+=1
                cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$(count) $gaps_type`
                push!(cmds,cmd)
            end
        end
        @showprogress "L=$(L) $exp_label" pmap(run,cmds)
        for i in 1:count
            # try
            open("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv","r") do input
                open("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv","a") do output
                    n = countlines(input)
                    record = 0
                    seekstart(input)
                    for k in 1:n
                        line = readline(input)
                        println(output,line)
                        record = record + 1
                    end
                    # println(record)
                end
            end
        #     catch
        #     end
        end
        for i in 1:count
            try
            rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
            catch
            end
        end
        for i in 1:count
            try
                open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv","r") do input
                    open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label).csv","a") do output
                        n = countlines(input)
                        record = 0
                        seekstart(input)
                        for k in 1:n
                            line = readline(input)
                            println(output,line)
                            record = record + 1
                        end
                        # println(record)
                    end
                end
            catch
            end
        end
        for i in 1:count
            try
                rm("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv")
            catch
            end
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
        # run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`) analysis of gaps is postponed for sake of time
    end
end
println("tasks completed.")

>>>>>>> cd1b2244bd3b8da94c96088688359f38b8af2f05
