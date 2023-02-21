# assume that the intensity data is stored in /intensity_exp
# "wt_15mM_salt" fold = 0 => /intensity_exp/intensity_wt0.csv
# "wt_150mM_salt" fold = 10 => /intensity_exp/intensity_lwt10.csv
# now we create a summary dataframe containing the relative intensity of each condition
# and the exp_label ∈ ["wt_15mM_salt", "wt_150mM_salt"]
using DataFrames, CSV
df = DataFrame(
    relative_intensity = Float64[],
    exp_label = String[],
    fold = Float64[],
    time = Float64[]
)

exp_dict = Dict(
    "wt_15mM_salt" => "./intensity_exp/intensity_wt",
    "wt_150mM_salt" => "./intensity_exp/intensity_lwt"
)

function exp_label_fold2filename(exp_label, fold)
    return exp_dict[exp_label] * string(fold) * ".csv"
end

folds = Int[0, 1, 4, 10, 25]
exp_labels = ["wt_15mM_salt", "wt_150mM_salt"]

for exp_label ∈ exp_labels, fold ∈ folds
    filename = exp_label_fold2filename(exp_label, fold)
    intensity = CSV.read(filename, DataFrame)
    relative_intensity = intensity.Normalized_intensity_increment .+ 1
    time = [5i for i in 1:length(relative_intensity)]
    data_df= DataFrame(
        relative_intensity = relative_intensity, 
        exp_label = repeat([exp_label], length(relative_intensity)), 
        fold = repeat([fold], length(relative_intensity)),
        time = time
        )
    df = vcat(df, data_df)
end

CSV.write("./intensity_exp/summary_intensity.csv", df)