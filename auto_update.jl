print("initializing...\n")
using Distributed
using DataFrames
using CSV

addprocs(32) # determined by the number of processors (cores)
@everywhere begin
    using ProgressMeter
end
@info "multi-threading initialized."
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$count.csv"
end


L = 5000

simupath = "$(pwd())/data_simu"
# initialize
if length(ARGS) == 1
    exp_label = ARGS[1]
else
    exp_label="wt_15mM_salt"
end
para_df = CSV.read("figs/sources/ini_para.csv",DataFrame)

init_label="init"
simu_label = "update"
gaps_type = "none"
it = 0
paras_names = ["k_on","k_off","v_open","v_close"]
figpath = "./figs"
requested_df=para_df[[para_df.exp_label[i] == exp_label for i in 1:length(para_df.L)],:]
k_on,k_off,k_open,k_close,α,β,L,exp_label = requested_df[1,:]

it_max = 15
# force resetting k_off
# k_off = 1e-6
current_loss = NaN
p₀=[k_on,k_off,k_open,k_close,current_loss]
# for para in paras_names
#     @show landscape=CSV.read("$figpath/sources/landscape_$(para)_$(exp_label)_$(init_label).csv",DataFrame)
#     p,i=findmin(landscape.error)
#     push!(p₀,landscape.para[i])
# end
P = []
para_df = DataFrame(k_on=Float64[],k_off=Float64[],v_open=Float64[],v_close=Float64[], loss=Float64[])
push!(para_df,p₀)
push!(P,p₀)
function simu_paras(p₀) # update the simulation parameters
    k_on,k_off,v_open,v_close,loss = p₀
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
    global k_ons = unique([modulate(k_on*(10.0^(i/8))) for i in -1:1:1])
    global k_offs = unique([modulate(k_off*(10.0^(i/8))) for i in -1:1:1])
    global v_opens = unique([modulate(v_open*(10.0^(i/8))) for i in -1:1:1])
    global v_closes = unique([modulate(v_close*(10.0^(i/8))) for i in -1:1:1])
end

function simu_paras(p₀,it) # update the simulation parameters
    k_on,k_off,v_open,v_close,loss = p₀
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
        global k_ons = unique([modulate(k_on*(10.0^(i/4))) for i in -1:1:1])
        global k_offs = unique([modulate(k_off*(10.0^(i/4))) for i in -1:1:1])
        global v_opens = unique([modulate(v_open*(10.0^(i/4))) for i in -1:1:1])
        global v_closes = unique([modulate(v_close*(10.0^(i/4))) for i in -1:1:1])
    end
end

@info "parameters initialized."

function get_paras(it)
    k_on, k_off, v_open, v_close, α, β, loss = CSV.read("figs/sources/result_$(exp_label)_$(simu_label)_$(it).csv",DataFrame)[1,:]
    return [k_on,k_off,v_open,v_close,loss]
end
## try to resume from previous fit
function get_current_it()
    it = 15
    while !isfile("$figpath/sources/plot_df_$(simu_label)_$(it)_$(exp_label).csv") && it > 0
        it -=1
    end
    return it
end
it = get_current_it()
if it > 0
    @info "resuming from iteration $(it)"
else 
    @info "starting from scratch"
    it = 0
end
# iteration
while it < it_max
    # resume previous fitting
    if it > 0 && length(P) < it 
        for j ∈ 1:it-1
            pₙ = get_paras(it)
            push!(P,pₙ)
            push!(para_df,pₙ)
        end
    end
    if it > 0
        pₙ = get_paras(it)
        if pₙ == P[it] # P[it] is the previous fitting result
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
    N = 100
    cmds = Array{Cmd,1}()
    count = 0
    @info "k_ons" (k_ons)
    @info "k_offs" (k_offs)
    @info "v_opens" (v_opens)
    @info "v_closes" (v_closes)
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
    run(`julia fit_base.jl $exp_label $(simu_label)_$(it) $(T1+T2)`)
    @info "updated parameters in iteration $it" para_df
end

## final stage
it = get_current_it()
p₁ = get_paras(it)
function final_paras(p₀) # update the simulation parameters
    k_on,k_off,v_open,v_close,loss = p₀
    global k_ons = unique([k_on for i in 1:1:1])
    global k_offs = unique([k_off for i in 1:1:1])
    global v_opens = unique([v_open for i in 1:1:1])
    global v_closes = unique([v_close for i in 1:1:1])
end
final_paras(p₁)
simu_label="fitted"
folds = [0,1,4,10,25,50]
Ls = [L]
T1 = 1800.0
T2 = 600.0
N = 50
n = ceil(Int,N/10)
gaps_type = "none"
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
    progress_pmap(run,cmds; 
        progress=Progress(length(cmds),desc="running for fitted parameters", 
                showspeed=true, color=:white)
                )
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
    # for i in 1:count
    #     try
    #         open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv","r") do input
    #             open("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label).csv","a") do output
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
    #     catch
    #     end
    # end
    # for i in 1:count
    #     try
    #         rm("$simupath/rsa_$(gaps_type)_gaps_$(exp_label)_$(simu_label)_$i.csv")
    #     catch
    #     end
    # end
    run(`julia evaluate.jl $exp_label $(simu_label)`)
    # run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)`)

