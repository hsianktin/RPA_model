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


L = 5000

simupath = "$(pwd())/data_simu"
# initialize
if length(ARGS) == 1
    exp_label = ARGS[1]
else
    exp_label="wt_150mM_salt"
end
para_df = CSV.read("figs/sources/para.csv",DataFrame)

init_label="init"
simu_label = "update"
gaps_type = "none"
it = 0
paras_names = ["k_on","k_off","v_open","v_close"]
figpath = "./figs"
requested_df=para_df[[para_df.L[i] == L && para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
k_on,k_off,k_open,k_close,α,β,L,exp_label = requested_df[1,:]

# force resetting k_off
# k_off = 1e-6

p₀=[k_on,k_off,k_open,k_close]
# for para in paras_names
#     @show landscape=CSV.read("$figpath/sources/landscape_$(para)_$(exp_label)_$(init_label).csv",DataFrame)
#     p,i=findmin(landscape.error)
#     push!(p₀,landscape.para[i])
# end
P = []
para_df = DataFrame(k_on=Float64[],k_off=Float64[],v_open=Float64[],v_close=Float64[])
push!(para_df,p₀)
push!(P,p₀)
function simu_paras(p₀) # update the simulation parameters
    k_on,k_off,v_open,v_close = p₀
    lower_limit = 1.0e-6
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
    global k_ons = unique([modulate(k_on*(10.0^(i/2))) for i in -1:1:1])
    global k_offs = unique([modulate(k_off*(10.0^(i/2))) for i in -1:1:1])
    global v_opens = unique([modulate(v_open*(10.0^(i/2))) for i in -1:1:1])
    global v_closes = unique([modulate(v_close*(10.0^(i/2))) for i in -1:1:1])
end

function simu_paras(p₀,it) # update the simulation parameters
    k_on,k_off,v_open,v_close = p₀
    lower_limit = 1.0e-6
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
    if it < 6
        global k_ons = unique([modulate(k_on*(10.0^(i))) for i in -1:1:1])
        global k_offs = unique([modulate(k_off*(10.0^(i))) for i in -1:1:1])
        global v_opens = unique([modulate(v_open*(10.0^(i))) for i in -1:1:1])
        global v_closes = unique([modulate(v_close*(10.0^(i))) for i in -1:1:1])
    else 
        global k_ons = unique([modulate(k_on*(10.0^(i/2))) for i in -1:1:1])
        global k_offs = unique([modulate(k_off*(10.0^(i/2))) for i in -1:1:1])
        global v_opens = unique([modulate(v_open*(10.0^(i/2))) for i in -1:1:1])
        global v_closes = unique([modulate(v_close*(10.0^(i/2))) for i in -1:1:1])
    end
end

print("parameters initialized.\n")
# iteration
while it < 15
    if it > 0
        pₙ = []
        for para in paras_names
            landscape=CSV.read("$figpath/sources/landscape_$(para)_$(exp_label)_$(simu_label)_$(it).csv",DataFrame)
            p,i=findmin(landscape.error)
            push!(pₙ,landscape.para[i])
        end
        if pₙ == P[it]
            println("optimization completed.")
            break
        else
            push!(para_df,pₙ)
            push!(P,pₙ)
        end
    else
        pₙ = p₀
    end
    global it += 1
    simu_paras(pₙ,it)
    folds = [0,1,4,10,25]
    Ls = [L]
    T1 = 1800.0
    T2 = 600.0
    N = 50
    cmds = Array{Cmd,1}()
    count = 0
    println(k_ons)
    println(k_offs)
    println(v_opens)
    println(v_closes)
    # while count < 200
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
    @showprogress 1 "running for iteration $(it) " pmap(run,cmds)
    try # clear old data
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it).csv")
        println("old data deleted")
    catch
        println("no old data detected")
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
    run(`julia evaluate.jl $exp_label $(simu_label)_$(it)`)
    @show para_df
end

## final stage
it = 15
while !isfile("$figpath/sources/landscape_$(paras_names[1])_$(exp_label)_$(simu_label)_$(it).csv")
    global it -=1
end
p₁ = []
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
simu_label="fitted"
folds = [0,1,4,10,25,50]
Ls = [L]
T1 = 1800.0
T2 = 600.0
N = 100
n = ceil(Int,N/10)
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
                                global count += 1
                                if !isfile(simu_fname(exp_label,simu_label,count))
                                    cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $n $(exp_label) $(simu_label)_$count $gaps_type`
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
    try
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
        println("previous results cleared")
    catch
        println("no previous results")
    end
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
    run(`julia evaluate.jl $exp_label $(simu_label)`)
    run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`)

