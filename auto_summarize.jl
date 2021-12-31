print("initializing...\n")
using DataFrames
using CSV
using ProgressMeter
function simu_fname(exp_label,simu_label,count)
    return "$simupath/rsa_plot_$(exp_label)_$(simu_label)_$(it)_$count.csv"
end


L = 5000

simupath = "$(pwd())/data_simu"

exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
para_df = DataFrame(k_on=Float64[],k_off=Float64[],v_open=Float64[],v_close=Float64[],α= Float64[],β=Float64[],L=Int64[],exp_label=String[])
paras_names = ["k_on","k_off","v_open","v_close","α","β"]
figpath = "./figs"
init_label="fitted"
simu_label = "update"
gaps_type = "none"

for exp_label in exp_labels
    it = 15
    while !isfile("$figpath/landscape/landscape_$(paras_names[1])_$(exp_label)_$(simu_label)_$(it).csv")
        it -= 1
    end
    # summarize all update data
    for i in 1:it
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
    for i in 1:it
        try
        rm("$simupath/rsa_plot_$(exp_label)_$(simu_label)_$i.csv")
        catch
        end
    end
    # run(`julia evaluate.jl $exp_label $(simu_label)`)
    run(`julia preprocess_data.jl $exp_label $(simu_label)`)
    run(`julia evaluate_over_preprocessed.jl $exp_label $(simu_label)`)
end

for exp_label in exp_labels
    p₀ = []
    for para in paras_names
        @show landscape=CSV.read("$figpath/landscape/landscape_$(para)_$(exp_label)_$(simu_label).csv",DataFrame)
        p,i=findmin(landscape.error)
        push!(p₀,landscape.para[i])
    end
    push!(para_df,vcat(p₀,[L,exp_label]))
end

CSV.write("$figpath/para_$(simu_label).csv",para_df)