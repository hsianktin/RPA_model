"""
This script is used to quickly summarize the sensitivity analysis results,
based on diff_exp_label_simu_label.csv generated by analysis script
"""


using CSV, DataFrames, Statistics

perturb_ids = [64, 107, 177, 202, 251, 316, 415, 489, 758, 794, 935, 940]
exp_labels = ["wt_15mM_salt", "wt_150mM_salt"]

# the standard paras for each experiment
function paras₀(exp_label)
    para_df = CSV.read("figs/para_fitted.csv", DataFrame)
    para_df = filter(row -> row.exp_label == exp_label, para_df)
    return [para_df.k_on[1], para_df.k_off[1], para_df.v_open[1], para_df.v_close[1], para_df.α[1], para_df.β[1]]
end


for exp_label in exp_labels
    perturbed_dfs = []
    for id in perturb_ids
        simu_label = "perturb_$id"
        df = CSV.read("data_simu/diff_$(exp_label)_$(simu_label).csv", DataFrame)
        push!(perturbed_dfs, df)
    end
    perturbed_df = vcat(perturbed_dfs...)
    k_on,k_off,v_open,v_close,α, β = paras₀(exp_label)

    perturbed_df_α = filter(row -> row.k_on == k_on && row.k_off == k_off && row.v_open == v_open && row.v_close == v_close && row.β == β, perturbed_df)
    # there are multiple rows with the same k_on, k_off, v_open, v_close, β
    # take the mean and std of the error
    perturbed_df_α = combine(groupby(perturbed_df_α, [:k_on, :k_off, :v_open, :v_close, :α, :β]), :diff => mean, :diff => std)
    CSV.write("figs/sources/landscape_α_$(exp_label)_perturb.csv", perturbed_df_α)

    # same for β
    perturbed_df_β = filter(row -> row.k_on == k_on && row.k_off == k_off && row.v_open == v_open && row.v_close == v_close && row.α == α, perturbed_df)
    perturbed_df_β = combine(groupby(perturbed_df_β, [:k_on, :k_off, :v_open, :v_close, :α, :β]), :diff => mean, :diff => std)
    CSV.write("figs/sources/landscape_β_$(exp_label)_perturb.csv", perturbed_df_β)
end

