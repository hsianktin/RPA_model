# based on existing parameters, obtain the relation between folds and covered DNA
print("initializing...\r")
using DataFrames
using Distributed
using CSV
addprocs(12)
@everywhere begin
    using ProgressMeter
end

para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# below is the old parameter obtained with L = 1000
## push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
## push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])

N = 100
T1 = 1800.0
T2 = 600.0
gaps_type = "none" # only obtain the coverage.
# simu_folds = [i for i in 0:2:100]
simu_folds = [0,1,4,10,25]
simu_label = "coverage_analysis" # only obtain the coverage
simupath = "$(pwd())/data_simu"

exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
print("initialized.....\n")

for exp_label in exp_labels
    L = 5000
    print("processing $exp_label...\n")
    requested_df=para_df[[para_df.L[i] == L && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
    k_on,k_off,k_open,k_close,α,β,L,exp_label,loss = requested_df[1,:]
    cmds = Array{Cmd,1}()
    count = 0
    for fold in simu_folds
        count += 1
        cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type`
        push!(cmds,cmd)
    end
    @showprogress "running simulation for $exp_label..." pmap(run,cmds)
    for i in 1:count
        
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
        
    end
    for i in 1:count
        try
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
        catch
        end
    end

    print("simulation for $exp_label completed.\n")
    include("evaluate_base.jl")
    initialize(exp_label,simu_label)
    print("\n$exp_label data loaded...\n")
    occupancy_plot(k_on,k_off,k_open,k_close,simu_folds,exp_label,"$(simu_label)_$L")
    print("$exp_label plot completed.\n")
    rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
end

# for exp_label in exp_labels
#     L = 1000
#     print("processing $exp_label...\n")
#     requested_df=para_df[[para_df.L[i] == L && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
#     k_on,k_off,k_open,k_close,α,β,L,exp_label = requested_df[1,:]
#     cmds = Array{Cmd,1}()
#     count = 0
#     for fold in simu_folds
#         count += 1
#         cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type`
#         push!(cmds,cmd)
#     end
#     @showprogress "running simulation for $exp_label..." pmap(run,cmds)
#     for i in 1:count
#         open("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv","r") do input
#             open("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv","a") do output
#                 n = countlines(input)
#                 record = 0
#                 seekstart(input)
#                 for k in 1:n
#                     line = readline(input)
#                     println(output,line)
#                     record = record + 1
#                 end
#                 # println(record)
#             end
#         end
#     end
#     for i in 1:count
#         try
#         rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
#         catch
#         end
#     end

#     print("simulation for $exp_label completed.\n")
#     include("evaluate_base.jl")
#     initialize(exp_label,simu_label)
#     print("\n$exp_label data loaded...\n")
#     occupancy_plot(k_on,k_off,k_open,k_close,simu_folds,exp_label,"$(simu_label)_$L")
#     print("$exp_label plot completed.\n")
#     rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")

# end