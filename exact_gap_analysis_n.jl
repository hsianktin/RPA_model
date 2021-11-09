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
exp_label = ""
simu_label = "wt_150mM_salt_unq"
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
Overall_Gaps = Array{Any,1}()
for _ in 1:5
    Gaps = Array{Float64,1}[]
    for i in 1:2401
        push!(Gaps,Array{Int16,1}())
    end
    push!(Overall_Gaps,Gaps)
end
function position_of_fold(fold)
    folds = [0,1,4,10,25]
    return findall(x-> x==fold,folds)[1]
end
function dbpath(exp_label,simu_label)
    path = "$simupath/rsa_exact_gaps_$(exp_label)_$(simu_label).csv"
    return path
end
function load(db)
    f = open(db,"r")
    n = countlines(f)  # length of files
    N = floor(Int,n/50) # number of samples
    seekstart(f)
    for i in 1:n
        # read simu_data.
        temp_line_1 = [parse(Float64,i) for i in split(readline(f),",")]
        global k_on,k_off,v_open,v_close,fold,L,T1,T2,_,t = temp_line_1[1:10]
        t = convert(Int,t)
        Overall_Gaps[position_of_fold(fold)][t] = [[convert(Int32,x) for x in temp_line_1[11:end]];Overall_Gaps[position_of_fold(fold)][t]]
        @printf("loading %s: %.3f\r",db,(i/n))
    end
    close(f)
    # data_frame = DataFrame(k_on=k_ons,k_off=k_offs,v_open=v_opens,v_close=v_closes,fold=data_folds,L=Ls,state_1=state_1s,state_2=state_2s)
end
load(dbpath(exp_label,simu_label))

T = [i for i in 0:1:2400]
T₀ = 1800
T₁ = 2400
using StatsPlots
dataframe = DataFrame(fold = Int[], size = Int[])
for fold in [0,1,4,10,25]
    data = Overall_Gaps[position_of_fold(fold)][T₁+1]
    for i in 1:length(data)
        push!(dataframe,[fold,data[i]])
    end
end

@df dataframe groupedhist(:size, group = :fold, bar_position = :dodge)
xlims!(0,40)
using CSV
CSV.write("./figs/exact_gaps_table_$(exp_label)_$(simu_label).csv",dataframe)
savefig("./figs/exact_gaps_hist_$(exp_label)_$(simu_label).png")
