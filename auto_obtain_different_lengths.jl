print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(12) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
print("multi-threading initialized.\n")
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$count.csv"
end

Lengths = [200,800,1000,5000]

simupath = "$(pwd())/data_simu"
# initialize
if length(ARGS) == 1
    exp_label = ARGS[1]
else
    exp_label="wt_150mM_salt"
end
init_label="init"
gaps_type = "exact"
it = 0
paras_names = ["k_on","k_off","v_open","v_close"]
figpath = "./figs"

function simu_paras(p₀) # update the simulation parameters
    k_on,k_off,v_open,v_close = p₀
    lower_limit = 1e-6
    upper_limit = 1e1
    function modulate(v) # fall between limits
        if v < lower_limit
            return lower_limit
        elseif v > upper_limit
            return upper_limit
        else
            return v
        end
    end
    global k_ons = unique([modulate(k_on*(10.0^i)) for i in -1:1:1])
    global k_offs = unique([modulate(k_off*(10.0^i)) for i in -1:1:1])
    global v_opens = unique([modulate(v_open*(10.0^i)) for i in -1:1:1])
    global v_closes = unique([modulate(v_close*(10.0^i)) for i in -1:1:1])
end
print("parameters initialized.\n")

## final stage
simu_label = "update"
if !isfile("$figpath/sources/landscape_$(paras_names[1])_$(exp_label)_$(simu_label)_$(it).csv") && it >= 0
    print("no fitting results found.\n Please run auto_summarize.jl\n")
    exit()
end

print("fitting parameters found")
p₁ = []
para_df = CSV.read("$figpath/para_$(simu_label).csv",DataFrame)
id = findfirst(para_df["exp_label"] == exp_label)
for para in paras_names
    @show landscape=CSV.read("$figpath/sources/landscape_$(para)_$(exp_label)_$(simu_label)_$(it).csv",DataFrame)
    p,i=findmin(landscape.error)
    push!(p₁,landscape.para[i])
end
function final_paras(p₀) # update the simulation parameters
    k_on,k_off,v_open,v_close = p₀
    global k_ons = unique([k_on for i in -1:1:1])
    global k_offs = unique([k_off for i in -1:1:1])
    global v_opens = unique([v_open for i in -1:1:1])
    global v_closes = unique([v_close for i in -1:1:1])
end
final_paras(p₁)
for L in Lengths
    simu_label="fitted_length_$L"
    folds = [0,1,4,10,25]
    Ls = [L]
    T1 = 1800.0
    T2 = 600.0
    N = 10
    gaps_type = "exact"
    cmds = Array{Cmd,1}()
    count = 0
    while count < 10*length(folds)
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
                                    cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type`
                                    push!(cmds,cmd)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # println("command length: $(length(cmds))")
    # use distributed storage for parallel execution.
    @showprogress 1 "running for fitted values " pmap(run,cmds)
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
    run(`julia evaluate_1.jl $exp_label $(simu_label)`)
    run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`)
end