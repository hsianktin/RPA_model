using DataFrames
using CSV
using Statistics
ids = [43, 133, 256, 653, 814, 881, 972]
exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
simu_label = "perturb"
data_path = "./figs/sources/"


parameters = ["k_on","k_off","v_open","v_close"]

function fname(parameter,simu_label,exp_label,id)
  return "$(data_path)landscape_$(parameter)_$(exp_label)_$(simu_label)_$(id).csv"
end

# run analysis
for exp_label in exp_labels
    for id in ids
        run(`julia evaluate_sensitivity.jl $exp_label $(simu_label)_$id`)
    end
end

# summarize data
for exp_label in exp_labels
    # exp_label = exp_labels[1] # test
    for parameter in parameters
        # parameter = parameters[1]
        para = Array{Float64,1}()
        para_flag = true
        errs = Array{Float64,1}[]
        for id in ids
            # id = ids[1]
            temp_df = CSV.read(fname(parameter,simu_label,exp_label,id),DataFrame)
            if para_flag
                for i in 1:length(temp_df.para)
                    push!(para,temp_df.para[i])
                end
                para_flag = false
            end
            push!(errs,temp_df.error)
        end
        μ_err = [mean([errs[j][i] for j in 1:length(errs)]) for i in 1:length(errs[1])]
        σ_err = [std([errs[j][i] for j in 1:length(errs)]) for i in 1:length(errs[1])]
        df = DataFrame(para=para,mean_err=μ_err,std_err=σ_err)
        CSV.write("$(data_path)/landscape_$(parameter)_$(exp_label)_$(simu_label).csv",df)
    end
end