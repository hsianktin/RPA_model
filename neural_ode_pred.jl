using CSV, DataFrames, Flux, Zygote, DifferentialEquations, Printf
df = CSV.read("data_simu/rsa_plot_landscape_mean_L_5000.csv", DataFrame)
# # for the sake of simplicity, we first consider 1D problem
df1 = df[[row for row in 1:nrow(df) if df.Id[row] == 0], :]
df = CSV.read("data_simu/rsa_plot_landscape_mean_L_5000.csv", DataFrame)

df2 = df[[row for row in 1:nrow(df) if df.Id[row] == 1], :]
para_df = df1[:, 1:5]
if df1[:, 1:5] != df2[:, 1:5]
    @warn "The parameters of the two dataframes are not the same"
end
device = cpu
# # we also select the data from t = 0 to t = 1800, 
# # to avoid time switching effect temporarily
Y_df = df1[:, 11:end-120]
X_df = df2[:, 11:end-120]

# Define the neural network (NN), which takes the current state and parameters as input
NN = Chain(Dense(7, 64, tanh), Dense(64, 2)) |> cpu
# set the weight matrix of NN to be 0 pointwise
θ,re = Flux.destructure(NN)

# read parameter from BSONN file 
using BSON
BSON.@load "trained_params.bson" θ

function dU╱dt(U, θ, t)
    # get rows of U as X, Kon, Koff, Kopen, Kclose, Fold
    X, Y, Kon, Koff, Kopen, Kclose, Fold = eachrow(U) 
    ∂ₜX_known = Kon .* (1 .- X) - (Koff .+ Kclose) .* X
    ∂ₜY_known = Koff .* X - Kopen .* Y
    NN_output = re(θ)(hcat(X, Y, log.(Kon), log.(Koff), log.(Kopen), log.(Kclose), Fold)|> transpose)
    ∂ₜX_NN, ∂ₜY_NN = eachrow(NN_output)
    # dX╱dt = Kon .* (1 .- X) .- (Koff .+ Kclose) .* X .+ ( (re(θ)(
    #     hcat(X, log.(Kon), log.(Koff), log.(Kopen), log.(Kclose), Fold)|> transpose  ) ) |> transpose) .* X .* α
    dX╱dt = ∂ₜX_known + ∂ₜX_NN .* X .* α
    dY╱dt = ∂ₜY_known + ∂ₜY_NN .* Y .* α 
    dP╱dt = zeros(Float32, (5, length(X))) |> cpu
    ∂ₜU = vcat(dX╱dt |> transpose, dY╱dt |> transpose , dP╱dt)
    return (∂ₜU |> Array |> cpu)
end
α = 1f-3 |> device

function batch_predict(U, θ)
    Prob = ODEProblem(dU╱dt, U, (0.0f0, 1800.0f0), θ)
    Sol = solve(Prob, Tsit5(), saveat=0:5:1800)
    # Sol.u # 361 Array, convert to (361, 100) Matrix 
    SolU = Sol.u
    return (
        hcat([u[1,:] for u in SolU]...) |> transpose |> cpu, 
        hcat([u[2,:] for u in SolU]...) |> transpose |> cpu)
end
function get_u0(row)
    x0 = X_df[row, 1]
    y0 = Y_df[row, 1]
    p0 = para_df[row, :] |> collect
    u0 = vcat(x0, y0 , p0)
    return u0
end

function get_xs(row)
    xs = X_df[row, :] |> collect
    return xs
end

function get_ys(row)
    ys = Y_df[row, :] |> collect
    return ys
end


train_data = [(get_u0(row), get_xs(row), get_ys(row)) for row in 1:nrow(X_df)]  
# convert data to Float32
train_data = [(Float32.(u0), Float32.(xs), Float32.(ys)) for (u0, xs, ys) in train_data]

U0 = hcat([u0 for (u0, xs, ys) in train_data]...) |> device
Xs = hcat([xs for (u0, xs, ys) in train_data]...) |> device
Ys = hcat([ys for (u0, xs, ys) in train_data]...) |> device

X̂s, Ŷs = batch_predict(U0, θ)

using PGFPlotsX
using LaTeXStrings
for row in 1:nrow(X_df)
    u0 = U0[:, row]
    xs = Xs[:, row]
    ys = Ys[:, row]
    x̂s = X̂s[:, row]
    ŷs = Ŷs[:, row]
    # # plot the prediction

    plt = @pgf Axis(
        {
            width = "3.4in",
            height= "3.4in",
            xlabel = L"t",
            ylabel = L"p(t)",
            # legend_pos = "north west",
            xmin = 0,
            xmax = 1800,
            ymin = 0,
            ymax = 1,
            title = "Neural ODE prediction",
        },
        Plot(
            {
                color = "red",
                no_marks = true,
                mark = "*",
                style = "dashed",
                },
            Table([0:5:1800, x̂s])
        ),
        LegendEntry(L"\hat{p}_{20}(t)"),
        Plot(
            {
                color = "blue",
                no_marks = true,
                mark = "*",
                style = "dashed",
                },
            Table([0:5:1800, ŷs])
        ),
        LegendEntry(L"\hat{p}_{30}(t)"),
        Plot(
            {
                color = "red",
                only_marks = true,
                mark_repeat = 10,
                mark = "o",
                },
            Table([0:5:1800, xs])
        ),
        LegendEntry(L"p_{20}(t)"),
        Plot(
            {
                color = "blue",
                only_marks = true,
                mark_repeat = 10,
                mark = "o",
                },
            Table([0:5:1800, ys])
        ),
        LegendEntry(L"p_{30}(t)"),
    )
    pgfsave("figs/neural_ode_samples/training_set_sample_$(row).pdf", plt)
end