using Distributed
using DataFrames
using CSV
addprocs(12) # determined by the number of processors (cores)
simupath = "$(pwd())/data_simu"

@everywhere begin
    using ProgressMeter
end
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_state_transition_$(exp_label)_$(simu_label)_$count.csv"
end

# 15mM wt
df = DataFrame(k_on = Float64[], k_off = Float64[], v_open = Float64[], v_close = Float64[])
simu_label = "general"
exp_label = "wt_150mM_salt"
k_ons = [10.0^(-3+k) for k in -0.25:0.25:0.25]
k_offs = [10.0^(-6+k) for k in -0.25:0.25:0.25]
v_opens = [10.0^(-3+k) for k in -0.25:0.25:0.25]
v_closes = [10.0^(-1.5+k) for k in -0.25:0.25:0.25]
folds = [0,1,4,10,25]
Ls = [1000]
T1 = 1800.0
T2 = 600.0
N = 20
cmds = Array{Cmd,1}()
count = 0
for k_on in k_ons
    for k_off in k_offs
        for v_open in v_opens
            for v_close in v_closes
                push!(df,(k_on,k_off,v_open,v_close))
                for L in Ls
                    for fold in folds
                        global count += 1
                        if !isfile(simu_fname(exp_label,simu_label,count))
                            cmd=`julia simu_state_transition_wt_15.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label)_$(simu_label)_$count`
                            push!(cmds,cmd)
                        end
                    end
                end
            end
        end
    end
end
CSV.write("./data_simu/history.csv",df)
println("command length: $(length(cmds))")
# use distributed storage for parallel execution.
@showprogress 1 "running simulation for $exp_label" pmap(run,cmds)
for i in 1:count
    try
    open("$simupath/rsa_state_transition_$(exp_label)_$(simu_label)_$i.csv","r") do input
        open("$simupath/rsa_state_transition_$(exp_label)_$(simu_label).csv","a") do output
            n = countlines(input)
            record = 0
            seekstart(input)
            for k in 1:n
                line = readline(input)
                println(output,line)
                record = record + 1
            end
            println(record)
        end
    end
    catch
    end
end
for i in 1:count
    try
    rm("$simupath/rsa_state_transition_$(exp_label)_$(simu_label)_$i.csv")
    catch
    end
end


# # 15mM wt
# df = DataFrame(k_on = Float64[], k_off = Float64[], v_open = Float64[], v_close = Float64[])
# simu_label = "crude"
# exp_label = "wt_150mM_salt"
# k_ons = [10.0^(k) for k in -6:2:-1]
# k_offs = [10.0^(k) for k in -6:2:-1]
# v_opens = [10.0^(k) for k in -4:2:0]
# v_closes = [10.0^(k) for k in -4:2:0]
# folds = [0,1,4,10,25]
# Ls = [1000]
# T1 = 1860.0
# T2 = 540.0
# N = 20
# cmds = Array{Cmd,1}()
# count = 0
# for k_on in k_ons
#     for k_off in k_offs
#         for v_open in v_opens
#             for v_close in v_closes
#                 push!(df,(k_on,k_off,v_open,v_close))
#                 for L in Ls
#                     for fold in folds
#                         global count += 1
#                         if !isfile(simu_fname(exp_label,simu_label,count))
#                             cmd=`julia simu_state_transition_wt_150.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label)_$(simu_label)_$count`
#                             push!(cmds,cmd)
#                         end
#                     end
#                 end
#             end
#         end
#     end
# end
# CSV.write("./data_simu/history.csv",df)
# println("command length: $(length(cmds))")
# # use distributed storage for parallel execution.
# @showprogress 1 "running simulation for $exp_label" pmap(run,cmds)
# for i in 1:count
#     try
#     open("$simupath/rsa_state_transition_$(exp_label)_$(simu_label)_$i.csv","r") do input
#         open("$simupath/rsa_state_transition_$(exp_label)_$(simu_label).csv","a") do output
#             n = countlines(input)
#             record = 0
#             seekstart(input)
#             for k in 1:n
#                 line = readline(input)
#                 println(output,line)
#                 record = record + 1
#             end
#             println(record)
#         end
#     end
#     catch
#     end
# end
# for i in 1:count
#     try
#     rm("$simupath/rsa_state_transition_$(exp_label)_$(simu_label)_$i.csv")
#     catch
#     end
# end

