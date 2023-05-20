# summarize the data from auto_update.jl
# create figs/para_fitted.csv
print("initializing...\n")
using DataFrames
using CSV
using DataFrames
using ProgressMeter
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$count.csv"
end


L = 5000

simupath = "$(pwd())/data_simu"

exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
para_df = DataFrame(k_on=Float64[],k_off=Float64[],v_open=Float64[],v_close=Float64[],α= Float64[],β=Float64[],L=Int64[],exp_label=String[], loss=Float64[])
paras_names = ["k_on","k_off","v_open","v_close","α","β"]
figpath = "./figs"
init_label="fitted"
simu_label = "fitted"
gaps_type = "none"

for exp_label in exp_labels
    run(`julia evaluate.jl $exp_label $(simu_label)`)
end

for exp_label in exp_labels
    # p₀ = []
    # for para in paras_names
    #     @show landscape=CSV.read("$figpath/sources/landscape_$(para)_$(exp_label)_$(simu_label).csv",DataFrame)
    #     p,i=findmin(landscape.error)
    #     push!(p₀,landscape.para[i])
    # end
    temp_df = CSV.read("figs/sources/loss_$(exp_label)_$(simu_label).csv",DataFrame)
    print(temp_df)
    k_on,k_off,v_open,v_close,α,β,loss = temp_df[1,:]
    push!(para_df,[k_on,k_off,v_open,v_close,α,β,L,exp_label,loss])
end

CSV.write("$figpath/para_$(simu_label).csv",para_df)