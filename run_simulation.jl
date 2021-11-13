using Distributed
using DataFrames
using CSV
addprocs(12) # determined by the number of processors (cores)
simupath = "$(pwd())/data_simu"

@everywhere begin
    using ProgressMeter
end
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$count.csv"
end

# 15mM wt
# LL = [300,800,1000,5000]
LL = [1000]
for L_used in LL
    df = DataFrame(k_on = Float64[], k_off = Float64[], v_open = Float64[], v_close = Float64[])
    simu_label = "update"
    exp_label = "wt_15mM_salt"
    gaps_type = "none"
    k_ons = [10.0^(-k) for k in 1:1:5]
    # K_D = [1,10,100,0.1,0.01]
    k_offs = [10.0^(-k) for k in 1:1:5]
    # k_offs = K_D.*1e1
    v_opens = [10.0^(-k) for k in 1:1:5]
    # v_opens = [10.0^(-i) for i in -1:0.5:1.5]
    v_closes = [10.0^(-k) for k in 1:1:5]

    folds = [0,1,4,10,25]
    Ls = [L_used]
    T1 = 1800.0
    T2 = 600.0
    N = 50
    cmds = Array{Cmd,1}()
    count = 0
    # while count < 200
        for k_on in k_ons
            for k_off in k_offs
                # k_off = k_on*K_D
                for v_open in v_opens
                    # for K_d in K_D
                        # v_close = K_d*v_open
                    for v_close in v_closes
                        push!(df,(k_on,k_off,v_open,v_close))
                        for L in Ls
                            for fold in folds
                                count += 1
                                if !isfile(simu_fname(exp_label,simu_label,count))
                                    cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type`
                                    push!(cmds,cmd)
                                end
                            end
                        end
                    end
                end
            end
        end
    # end
    CSV.write("./data_simu/history.csv",df)
    println("command length: $(length(cmds))")
    # use distributed storage for parallel execution.
    @showprogress 1 "running simulation for $exp_label" pmap(run,cmds)
    for i in 1:count
        try
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
        catch
        end
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
end

# # # 150mM wt
# # LL = [1000]
# for L_used in LL
#     df = DataFrame(k_on = Float64[], k_off = Float64[], v_open = Float64[], v_close = Float64[])
#     # simu_label = "length_$(L_used)"
#     simu_label = "final_length_$(L_used)"
#     exp_label = "wt_150mM_salt"
#     gaps_type = "cumulative"
#     k_ons = [10.0^(-5) for k in 1:1:5]
#     # K_D = [1,10,100,0.1,0.01]
#     k_offs = [10.0^(-3) for k in 1:1:2]
#     # k_offs = K_D.*1e1
#     v_opens = [10.0^(-2) for k in 1:1:1]
#     # v_opens = [10.0^(-i) for i in -1:0.5:1.5]
#     v_closes = [10.0^(4) for k in 1:1:1]

#     folds = [0,1,4,10,25]
#     Ls = [L_used]
#     T1 = 1800.0
#     T2 = 600.0
#     N = 10
#     cmds = Array{Cmd,1}()
#     count = 0
#     # while count < 200
#         for k_on in k_ons
#             for k_off in k_offs
#                 # k_off = k_on*K_D
#                 for v_open in v_opens
#                     # for K_d in K_D
#                         # v_close = K_d*v_open
#                     for v_close in v_closes
#                         push!(df,(k_on,k_off,v_open,v_close))
#                         for L in Ls
#                             for fold in folds
#                                 count += 1
#                                 if !isfile(simu_fname(exp_label,simu_label,count))
#                                     cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type`
#                                     push!(cmds,cmd)
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     # end
#     CSV.write("./data_simu/history.csv",df)
#     println("command length: $(length(cmds))")
#     # use distributed storage for parallel execution.
#     @showprogress 1 "running simulation for $exp_label" pmap(run,cmds)
#     for i in 1:count
#         try
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
#                 println(record)
#             end
#         end
#         catch
#         end
#     end
#     for i in 1:count
#         try
#         rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
#         catch
#         end
#     end
#     for i in 1:count
#         try
#             open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv","r") do input
#                 open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label).csv","a") do output
#                     n = countlines(input)
#                     record = 0
#                     seekstart(input)
#                     for k in 1:n
#                         line = readline(input)
#                         println(output,line)
#                         record = record + 1
#                     end
#                     println(record)
#                 end
#             end
#         catch
#         end
#     end
#     for i in 1:count
#         try
#             rm("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv")
#         catch
#         end
#     end
# end

run(`julia run_analysis.jl`)