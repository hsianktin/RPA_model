# find the local sensitivity with respect to diffusion
print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(6) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
print("multi-threading initialized.\n")

function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$count.csv"
end

para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
# push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])
# CSV.write("figs/sources/para.csv",para_df)
N = 5
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"
simupath = "$(pwd())/data_simu"
folds = [0,1,4,10,25,50]
simu_folds = [0,1,4,10,25,50]
Ds = [10.0^(i) for i in 3:3] # range of diffusion
# function pert_paras(p₀,index) # update the simulation parameters
#     k_on,k_off,v_open,v_close = p₀
#     lower_limit = 1e-6
#     upper_limit = 1e1
#     function modulate(v) # fall between limits
#         if v < lower_limit
#             return lower_limit
#         elseif v > upper_limit
#             return upper_limit
#         else
#             return v
#         end
#     end
#     if index == 1
#         global k_ons = unique([modulate(k_on*(10.0^(i))) for i in -2:0.5:2])
#     elseif index == 2
#         global k_offs = unique([modulate(k_off*(10.0^(i))) for i in -2:0.5:2])
#     elseif index == 3
#         global v_opens = unique([modulate(v_open*(10.0^(i))) for i in -2:0.5:2])
#     elseif index == 4
#         global v_closes = unique([modulate(v_close*(10.0^(i))) for i in -2:0.5:2])
#     end
# end
simupath="./data_simu"
L = 500
simu_label = "perturb_diffusion_$(rand([i for i in 1000:2000]))"
# simu_label = "perturb_diffusion_1573"
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
for exp_label in exp_labels
    cmds = Array{Cmd,1}()
    print("processing $exp_label...\n")
    # requested_df=para_df[[para_df.L[i] == L && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L )],:]
    requested_df=para_df[[para_df.exp_label[i] == exp_label for i in 1:length(para_df.L )],:]
    k_on,k_off,k_open,k_close,α,β,_,exp_label,loss = requested_df[1,:]
    count = 0
    global k_ons = [k_on]
    global k_offs = [k_off]
    global v_opens = [k_open]
    global v_closes = [k_close]
    # pert_paras([k_on,k_off,k_open,k_close],i)
    for k_on in k_ons
        for k_off in k_offs
            # k_off = k_on*K_D
            for v_open in v_opens
                # for K_d in K_D
                    # v_close = K_d*v_open
                for v_close in v_closes
                    for D in Ds
                    # push!(df,(k_on,k_off,v_open,v_close))
                        for fold in folds
                            count += 1
                            if !isfile(simu_fname(exp_label,simu_label,count))
                                cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type $D`
                                push!(cmds,cmd)
                            end
                        end
                    end
                end
            end
        end
    end
    @showprogress 1 "running for $exp_label " pmap(run,cmds)
    try # clear old data
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
        println("old data deleted")
    catch
        println("no old data detected")
    end
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
    run(`julia evaluate_diffusion.jl $exp_label $(simu_label)`)
    run(`julia exact_gap_analysis_diffusion_selective.jl $exp_label $simu_label`)
end
