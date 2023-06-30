using DataFrames, CSV
using Pipe
intensity_df = CSV.read(
    "./intensity_exp/summary_intensity.csv", 
    DataFrame
)

# access the experimental data
function access_intensity_trace(exp_label, fold)
    query_df = @pipe intensity_df |>
        filter(row -> row.exp_label == exp_label, _) |>
        filter(row -> row.fold == fold, _)
    relative_intensity = query_df.relative_intensity[1:481]
    time = query_df.time[1:481]
    return relative_intensity, time
end

# access the simulated data (modified in evaluate_base.jl)
# access_intensity_trace("wt_15mM_salt", 0)