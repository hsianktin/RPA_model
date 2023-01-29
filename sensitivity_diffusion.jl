# find the local sensitivity with respect to diffusion
@info "initializing..."
using Distributed
using DataFrames
using CSV
using Printf
using Random: shuffle
if nprocs() == 1
    addprocs([("rog-ld03", 24),("ryzen",16)]) # determined by the number of processors (cores)
    # addprocs(15)
end
@everywhere begin
    using ProgressMeter
end
@info "multi-threading initialized."

function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$count.csv"
end

para_df = CSV.read("figs/para_fitted.csv",DataFrame)
# push!(para_df,[1e-5,1e-3, 1e-4,1e-3,2,3.2,1000,"wt_15mM_salt"])
# push!(para_df,[1e-5,1e-2, 1e-4,1e-2,2,3.2,1000,"wt_150mM_salt"])
# CSV.write("figs/sources/ini_para.csv",para_df)
N = 100
dry = false
L = 5000

if length(ARGS) > 1
    if ARGS[1] == "--dry"
        # execute `julia sensitivity_diffusion.jl -- --dry` for testing
        N = 100
        dry = true
        L = 100
        @info "dry run mode for testing"
    end
end
T1 = 1800.0
T2 = 600.0
gaps_type = "none"
simupath = "$(pwd())/data_simu"
folds = [0,1,4,10,25,50]
simu_folds = [0,1,4,10,25,50]
Ds = [10.0^(i) for i in -5:1:0] # range of diffusion
A = [i/10 for i in 5:1:20]

nₘₐₓ = 5
simupath="./data_simu"

# set headless mode for plotting
ENV["GKSwstype"] = "100"


# add functionality to save the simulation label
# and to resume the simulation from the last saved label
if isfile("tmp/sensitivity_diffusion_simu_label.txt")
    @info "found existing simulation label, resuming..."
    simu_label = readlines("tmp/sensitivity_diffusion_simu_label.txt")[1]
    @info "simu_label = $simu_label"
else
    @info "no existing simulation label, generating a new one..."
    simu_label = "diffusion_$(rand([i for i in 1000:2000]))"
    @info "simu_label = $simu_label"
    open("tmp/sensitivity_diffusion_simu_label.txt","w") do f
        write(f,"$simu_label")
    end
end
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
for exp_label in exp_labels
    cmds = Array{Cmd,1}()
    @info "processing $exp_label..."
    if isfile("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
        @info "found existing data, skipping..."
        continue
    end
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
                            if N ≤ nₘₐₓ
                                if !isfile(simu_fname(exp_label,simu_label,count))
                                    cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label)_$count $gaps_type $D`
                                    push!(cmds,cmd)
                                end
                                count += 1
                            else
                                iters = Int(floor(N/nₘₐₓ))
                                for i in 1:iters
                                    if !isfile(simu_fname(exp_label,simu_label,count))
                                        cmd=`julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $nₘₐₓ $(exp_label) $(simu_label)_$count $gaps_type $D`
                                        push!(cmds,cmd)
                                    end
                                    count += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # @showprogress 1 "running simulation for $exp_label " pmap(run,cmds)
    progress_pmap(run, shuffle(cmds); # shuffle to get an accurate estimate of the time remaining
        progress=Progress(length(cmds), 
                barglyphs=BarGlyphs("[=> ]"), 
                showspeed=true, 
                desc="running simulation for $exp_label",
                color=:white,
                )
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
    include("evaluate_base_diffusion.jl")
    initialize(exp_label,simu_label)
    length_of_paras = length(index_p.k_on)
    @showprogress 1 "analyzing $exp_label..." for i in 1:length_of_paras
        local k_on,k_off,v_open,v_close,D1 = index_p[i,:]
        evaluate(df,k_on,k_off,v_open,v_close,D1,folds)
        # @printf("analyzing progress: %.2f\r",i/length_of_paras)
    end
    analyze(df,"D1",exp_label, simu_label)
    # analyze(df,"α")
    # analyze(df,"β")
    # CSV.write("./data_simu/analyze_$(exp_label)_$(simu_label).csv")
    minval,index_min = findmin(df.diff)
    @info "minimum loss = $minval"
    k_on,k_off,v_open,v_close,D1,α,β,d=df[index_min,:]
    @info @sprintf("minimum point: k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,D1=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,D1,α,β)
    # k_on,k_off,v_open,v_close,α,β= 1e-5,1e-4,1e-2,1e-3,1,2
    ensemble_plot(k_on,k_off,v_open,v_close,D1,folds,α,β; simu_label = simu_label, exp_label=exp_label)
    yaxis!(:flip)
    ylims!(-2.0,2.0)
    # xlims!(1500,2400)
    @info "saving ensemble_plot at ./figs/rsa_state_transition_$(exp_label)_$simu_label.svg"
    savefig("./figs/rsa_state_transition_$(exp_label)_$simu_label.svg")
    # microstates_plot(k_on,k_off,v_open,v_close,D1,folds, exp_label, simu_label)
    
    f = open("./figs/sources/loss_$(exp_label)_$(simu_label).csv","w")
    write(f,"$minval")
    close(f)
    # run(`julia exact_gap_analysis_diffusion_selective.jl $exp_label $simu_label`)
    if dry 
        # delete all generated files
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label).csv")
        rm("./figs/rsa_state_transition_$(exp_label)_$simu_label.svg")
        rm("./figs/sources/loss_$(exp_label)_$(simu_label).csv")
    end
end

# complete simulation, delete "tmp/sensitivity_diffusion_simu_label.txt" 
rm("./tmp/sensitivity_diffusion_simu_label.txt")

