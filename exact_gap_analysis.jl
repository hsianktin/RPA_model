using DataFrames: DataAPI
using Base: Float64
using Plots: push!, position
# Evaluating Gaps...
# Only suitable in the case where all data come from the same parameters
using Statistics
using CSV
using Plots
using DataFrames
using DelimitedFiles
using Printf
if length(ARGS) == 1
    exp_label = ""
    simu_label = ARGS[1]
    # simu_label = "wt_15mM_salt"
elseif length(ARGS) == 2
    exp_label=ARGS[1]
    println("exp_label=$exp_label")
    simu_label=ARGS[2]
    println("simu_label=$simu_label")
else
    exp_label = "wt_15mM_salt"
    simu_label = "length_1000"
end

simupath = "$(pwd())/data_simu"
mainpath = pwd()
exppath = "$mainpath/data_exp"
figpath = "$mainpath/figs"
k_ons = Array{Float64,1}()
k_offs = Array{Float64,1}()
v_opens = Array{Float64,1}()
v_closes = Array{Float64,1}()
data_folds = Array{Float64,1}() # used for recording data
Ls = Array{Float64,1}()
T1s = Array{Float64,1}()
T2s = Array{Float64,1}()
Ns = Array{Float64,1}()

function dbpath(exp_label,simu_label)
    path = "$simupath/rsa_exact_gaps_$(exp_label)_$(simu_label).csv"
    return path
end
function load(db)
    df = DataFrame(time = Int[], fold = Int[], size = Int[], id = Int[])
    f = open(db,"r")
    n = countlines(f)  # length of files
    N = floor(Int,n/50) # number of samples
    seekstart(f)
    id = 1
    for i in 1:n
        # read simu_data.
        temp_line_1 = [parse(Float64,i) for i in split(readline(f),",")]
        if i > 1 && temp_line_1[10] < t
            id += 1
        end
        global k_on,k_off,v_open,v_close,fold,L,T1,T2,_,t = temp_line_1[1:10]
        t = convert(Int,t)
        for s in [convert(Int32,x) for x in temp_line_1[11:end]]
            push!(df, [t, fold, s, id])
        end
        @printf("loading %s: %.3f\r",db,(i/n))
    end
    close(f)
    return df
    # data_frame = DataFrame(k_on=k_ons,k_off=k_offs,v_open=v_opens,v_close=v_closes,fold=data_folds,L=Ls,state_1=state_1s,state_2=state_2s)
end
df = load(dbpath(exp_label,simu_label))

T = [i for i in 0:1:2400]
T₀ = 1800
T₁ = 2400
using StatsPlots

@df df[df.time .== 2401,:] groupedhist(:size, group = :fold, bar_position = :dodge)
title!("Count=$(maximum(df.id)/length(unique(df.fold)))")
xlims!(0,40)
using CSV
CSV.write("./figs/exact_gaps_table_$(exp_label)_$(simu_label).csv",df[findall(x->x∈[30*60+1,31*60+1,32*60+1,40*60+1],df.time),:])
savefig("./figs/exact_gaps_hist_$(exp_label)_$(simu_label).png")
