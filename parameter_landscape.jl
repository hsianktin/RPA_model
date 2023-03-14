print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(8) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
@info "workers initialized."

function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$count.csv"
end


L = 500

simupath = "$(pwd())/data_simu"
# initialize
if length(ARGS) == 1
    exp_label = ARGS[1]
else
    exp_label="wt_15mM_salt"
end

simu_label = "landscape"
gaps_type = "none"
it = 0
paras_names = ["k_on","k_off","v_open","v_close"]
figpath = "./figs"
requested_df=para_df[[para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]

k_ons = [10.0^i for i in -6:0.5:0]
k_offs = [10.0^i for i in -6:0.5:0]
v_opens = [10.0^i for i in -6:0.5:0]
v_closes = [10.0^i for i in -6:0.5:0]
folds = [0,1,4,10,25]
Ls = [L]
T1 = 1800.0
T2 = 600.0
N = 100
cmds = Array{Cmd,1}()
count = 0

for k_on in k_ons
    for k_off in k_offs
        # k_off = k_on*K_D
        for v_open in v_opens
            # for K_d in K_D
                # v_close = K_d*v_open
            for v_close in v_closes
                # push!(df,(k_on,k_off,v_open,v_close))
                for L in Ls
                    for fold in folds
                        count += 1
                        if !isfile(simu_fname(exp_label,simu_label,count))
                            cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$(it)_$count $gaps_type`
                            push!(cmds,cmd)
                        end
                    end
                end
            end
        end
    end
end
# end
# println("command length: $(length(cmds))")
# use distributed storage for parallel execution.
progress_pmap(run,cmds; 
progress=Progress(length(cmds),desc="running for iteration $(it)", 
        showspeed=true, color=:white)
        )
try # clear old data
rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it).csv")
@info "old data deleted and removed"
catch
@info "no old data detected, saving..."
end
for i in 1:count
try
open("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$i.csv","r") do input
    open("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it).csv","a") do output
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
rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$i.csv")
catch
end
end
for i in 1:count
try
    open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$(it)_$i.csv","r") do input
        open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$(it).csv","a") do output
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
    rm("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$(it)_$i.csv")
catch
end
end
