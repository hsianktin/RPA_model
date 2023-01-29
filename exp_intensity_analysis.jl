using CSV,Pipe,DataFrames
using Statistics
using LinearAlgebra

intensity_path(exp_label,fold) = "intensity_exp/$(exp_label)$(fold)"

function load_intensity_df(exp_label,fold)
    path = intensity_path(exp_label,fold)
    df = DataFrame(
        time = Float64[],
        intensity = Float64[],
        fold = Int[],
        exp_label = String[],
        id = Int[],
    )
    for (i, file) in enumerate(readdir(path))
        file_path = joinpath(path,file)
        print(i)
        println(file_path)
        file_df = CSV.read(file_path, DataFrame; header=["intensity"])
        for j in 1:nrow(file_df)
            t = (j-0)*5
            intensity = file_df.intensity[j]
            push!(df, [t, intensity, fold, exp_label, i])
        end
    end
    return df
end


exp_label="exp_data_wt"
fold = 0
intensity_df = load_intensity_df(exp_label,fold)

time_ref = [355,356,357,477,478,479].*5.0

# linear regression
# summarizing function
# √(V::Vector) = sqrt.(V) 
∑(x) = sum(x) 
function lm(X,Y)
    n = length(Y)
    Ȳ = ∑(Y)/n
    if length(size(X))>1
        p = size(X)[2]
    else
        p = 1
    end
    β̂ = (X'*X)^(-1)*X'*Y
    σ̂ = norm(Y - X*β̂)/√(n-p)
    σ̂² = σ̂^2
    dof_residual = n-p
    covβ̂ = σ̂²*(X'*X)^(-1)
    TSS = ∑((Y.-Ȳ).^2) # total sum of squares
    SSE = ∑((Y-X*β̂).^2) # sum of squares of errors
    R² = 1-SSE/TSS
    R²_adj = 1-(SSE/dof_residual)/(TSS/(n-1))
    return LM(X,Y,n,β̂,σ̂²,p,covβ̂,TSS,SSE,R²,R²_adj)
end

struct LM
    # this struct include all the information about the linear model
    X
    Y
    n
    β̂
    σ̂²
    p
    covβ̂
    TSS
    SSE
    R²
    R²_adj
end

function lm(df::DataFrame)
    n = nrow(df)
    X = ones(n,2)
    X[:,2] = df.time
    Y = df.intensity
    return lm(X,Y)
end

df_wt = @pipe load_intensity_df("exp_data_wt",fold) |>
filter(row -> row.time ∈ time_ref, _)
@pipe df_wt |> filter(row -> row.time == 2390.0, _) |> std(_.intensity)

lm_wt = @pipe load_intensity_df("exp_data_wt",fold) |>
    filter(row -> row.time ∈ time_ref, _) |>
    lm 
PrT(lm_wt)
PrL2400(lm_wt)
lm_lwt = @pipe load_intensity_df("exp_data_lwt",fold) |>
    filter(row -> row.time ∈ time_ref, _) |> lm
PrT(lm_lwt)
PrL2400(lm_lwt)
using Distributions
function PrT(linearmodel)
    # T test
    return 2(1-cdf(TDist(linearmodel.n-linearmodel.p),
    abs(linearmodel.β̂[2]/√(linearmodel.covβ̂[2,2]))))
end

function PrL2400(linearmodel)
    v = [1,2400]
    ŷ = ∑(v.*linearmodel.β̂)
    println("ŷ = $(ŷ)")
    σ̂ = √(v'*linearmodel.covβ̂*v)
    println("σ̂ = $(σ̂)")
    return 2(1-cdf(TDist(linearmodel.n-linearmodel.p),abs(ŷ/σ̂)))
end


lm_1.β̂[2]/√(lm_1.covβ̂[2,2])