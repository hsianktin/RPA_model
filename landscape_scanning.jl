@info "initializing..."
using Distributed
using DataFrames
using CSV

addprocs(12) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
@info "workers initialized."

function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$count.csv"
end


L = 50000

simupath = "$(pwd())/data_simu"
# initialize
if length(ARGS) == 1
    exp_label = ARGS[1]
else
    exp_label=""
end

simu_label = "landscape_L_$(L)"
gaps_type = "none"
figpath = "./figs"

k_ons = [10.0^i for i in -3:1:-1]
k_offs = [10.0^i for i in -6:1:-3]
v_opens = [10.0^i for i in -6:1:-1]
v_closes = [10.0^i for i in -6:1:-1]
folds = [1]
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
                        global count += 1
                        if !isfile(simu_fname(exp_label,simu_label,count))
                            cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type 0.0 5`
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
# randomize cmds for better speed estiamtes
using Random
cmds = shuffle(cmds)
progress_pmap(run,cmds; 
progress=Progress(length(cmds),desc="running for the landscape", 
        showspeed=true, color=:white)
        )
try # clear old data
rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
@info "old data deleted and removed"
catch
@info "no old data detected, saving..."
end

# first print header to rsa_plot_$(exp_label)_$(simu_label).csv
cols = [
    "k_on",
    "k_off",
    "v_open",
    "v_close",
    "fold",
    "L",
    "T1",
    "T2",
    "D1",
    "Id",
]

dt = 5
ts = [string(i) for i in 0:dt:2400] 
cols = [cols; ts]
# concatenate the columns by ","
header = join(cols, ",")
# write the header to the file
open("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv","a") do output
    println(output,header)
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

## Read the CSV file and recompress it, delete the original file
using DataFrames, CSV
df = CSV.read("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv", DataFrame)
CSV.write("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv.gz",df; compress=true)
rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")


# when completed, send a notification
# to my email address
# using ssmtp -t 
# to send email
email_addr = "xiangting.li@ucla.edu"
email_body = """
From: Automata <$(email_addr)>
To: Xiangting <$(email_addr)>
Subject: RSA landscape simulation completed

A total of  $(count) simulations have been completed.
The results are saved in $(simupath)/rsa_plot_$(exp_label)_$(simu_label).csv

Automata
"""

# save email body to a file
open("email_body_$(exp_label)_$(simu_label).txt","w") do f
    println(f,email_body)
end
# send email
run(pipeline(`ssmtp -t`, stdin = "email_body_$(exp_label)_$(simu_label).txt"))
# remove email body file
rm("email_body_$(exp_label)_$(simu_label).txt")