# train a simple neural ODE model for the dynamics of mean occupancy.
using CSV 
using DataFrames
using Flux 
using DiffEqFlux,Zygote
using DifferentialEquations, Optimization, OptimizationOptimJL
using CUDA
using Printf
# set up device by checking the CUDA version
# if CUDA.functional()
#     device = gpu
# else
#     device = cpu
# end
device = cpu
CUDA.allowscalar(false)
# set allow scalar indexing to be false 
# to avoid the error of scalar indexing on GPU


# load data 
df = CSV.read("data_simu/rsa_plot_landscape_mean_L_5000.csv", DataFrame)
# # for the sake of simplicity, we first consider 1D problem
df1 = df[[row for row in 1:nrow(df) if df.Id[row] == 0], :]
df = CSV.read("data_simu/rsa_plot_landscape_mean_L_5000.csv", DataFrame)

df2 = df[[row for row in 1:nrow(df) if df.Id[row] == 1], :]
para_df = df1[:, 1:5]
if df1[:, 1:5] != df2[:, 1:5]
    @warn "The parameters of the two dataframes are not the same"
end
# # we also select the data from t = 0 to t = 1800, 
# # to avoid time switching effect temporarily
Y_df = df1[:, 11:end-120]
X_df = df2[:, 11:end-120]
# design the neural ODE model
# let x(t) denote the data, p denote the parameters 
# and f(x(t), p) denote the dynamics of x(t)
# dx/dt = f(x(t), p)
# we consider some known information,
#  f(x,p) = p[1] * (1 -x ) - (p[2] + p[4])x + NN(θ, p, x)
# for each row of X_df, we have 
# ts = [0:5:1800] |> collect 
# xs = X_df[row, :] |> collect

# Define the neural network (NN), which takes the current state and parameters as input
NN = Chain(Dense(7, 64, tanh), Dense(64, 2)) |> cpu
# set the weight matrix of NN to be 0 pointwise
θ,re = Flux.destructure(NN)
# Set the scale of the initial output of the NN 
# (this is optional but can help with training)

α = 1f-2 |> device

# Define the ODE right hand side function including the known dynamics and the NN
# this formulation is bad, because it will not be able to do scalar indexing on GPU
function du╱dt(u,θ,t) # θ is the parameters of the NN
    x, kon, koff, kopen, kclose, fold = u
    dx╱dt = kon * (1 - x) - (koff + kclose) * x +  cpu(re(θ)(
        vcat(x, log.([kon,koff,kopen,kclose,fold])) |> device)[1])  * x * α 
    dx╱dt=  Float32(dx╱dt)
    dp╱dt = zeros(Float32,5)
    du╱dt = vcat(dx╱dt, dp╱dt)
    return du╱dt 
end


# Wrap the ODE function into a Neural ODE
# neuralode = NeuralODE(du╱dt, (0.0, 1800.0), Tsit5(), saveat=0:5:1800)
# get input to neural ODE

# du╱dt(u,θ, 0.0) # θ is the parameters of the NN
# x, kon, koff, kopen, kclose, fold = u
# output = (vcat(x, log.([kon,koff,kopen,kclose,fold]))) |> device |> re(θ) 

# du╱dt(u0, θ, 0.0)
function predict(u0)
    node_prob =  ODEProblem(du╱dt, u0, (0.0f0, 1800.0f0), θ)
    us = solve(node_prob, Tsit5() , saveat=0:5:1800).u
    # implicit methods doesn't support backpropagation
    xs = [u[1] for u in us]
    return xs
end

# train data = [(get_u0(row), get_xs(row)) for row in 1:nrow(X_df)]
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
# Define the loss function as the sum of squared errors between the true and predicted values
loss(u0, xs, ys) = sum(abs, predict(u0) .- [xs;ys]) / (2length(xs))

U0 = hcat([u0 for (u0, xs, ys) in train_data]...) |> device
Xs = hcat([xs for (u0, xs, ys) in train_data]...) |> device
Ys = hcat([ys for (u0, xs, ys) in train_data]...) |> device

U = U0
# parallel computing of du╱dt()
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

U0 = hcat([u0 for (u0, xs) in train_data]...) |> cpu # 6 × N_sample Matrix
Xs = hcat([xs for (u0, xs) in train_data]...) |> cpu # 361 × N_sample Matrix

# U0 = rand(Float32, (6, 2)) |> cpu
# Xs = rand(Float32, (361, 2)) |> cpu


# dU╱dt(U, θ, 0.0) |> typeof # LinearAlgebra.Transpose{Float64, cpu{Float64, 2, CUDA.Mem.DeviceBuffer}}
# typeof(U)
# U + dU╱dt(U, θ, 0.0) .* 0.01

function batch_predict(U, θ)
    Prob = ODEProblem(dU╱dt, U, (0.0f0, 1800.0f0), θ)
    Sol = solve(Prob, Tsit5(), saveat=0:5:1800)
    # Sol.u # 361 Array, convert to (361, 100) Matrix 
    SolU = Sol.u
    return (
        hcat([u[1,:] for u in SolU]...) |> transpose |> cpu, 
        hcat([u[2,:] for u in SolU]...) |> transpose |> cpu)
end



batch_predict(U0, θ) |> typeof # 361 × N_sample Matrix

function batch_loss(U0, θ, Xs, Ys) 
    X̂s, Ŷs = batch_predict(U0, θ)
    return sum(abs, X̂s .- Xs) / (2length(Xs)) + sum(abs, Ŷs .- Ys) / (2length(Ys))
end
α = 0.0f0 |> cpu
@time batch_loss(U0, θ, Xs, Ys)

# test gradient of loss function
@time ∇θ = gradient(Flux.params(θ)) do
    batch_loss(U0, θ, Xs, Ys)
end

α = 0 |> device
@info "mean-field ODE loss" batch_loss(U0, θ, Xs, Ys)

α = 1f-3 |> device
@info "untrained neural ODE loss with α = $α" batch_loss(U0, θ, Xs, Ys)

# train 
opt = ADAMW(0.01)

# rearrange train_data by random permutation
# load packages for shuffle

using ProgressBars
n_epochs = 10000
p = ProgressBar(1:n_epochs)
losses = []
for iter in p
    # iter = 1
    ∇θ = gradient(Flux.params(θ)) do
        batch_loss(U0, θ, Xs, Ys)
    end
    current_loss = batch_loss(U0, θ, Xs, Ys)
    Flux.update!(opt, Flux.params(θ), ∇θ)
    set_multiline_postfix(p, "loss = $current_loss")
    push!(losses, current_loss)
end
    
@info "trained neural ODE loss with α = $α" batch_loss(U0, θ, Xs, Ys)

# save trained parameters and restructure function 
using BSON 
BSON.@save "trained_params.bson" θ
# save loss
CSV.write("loss.csv", DataFrame(loss=losses), writeheader=false)

# use Term plot to print loss vs. iteration
using UnicodePlots
using Term
panel(plot; kw...) = Panel(string(plot, color=true); fit=true, kw...)

print(
    panel(lineplot(Array{Float32}((losses));); title="NN losses")
)


# loss()

# # compute the gradient wrt θ
# ∇θ = gradient(Flux.params(θ)) do
#     @show loss()
# end

# θ -= ∇θ[θ]
# loss()
# θ
# Flux.update!(opt, Flux.params(θ), ∇θ)
# θ
# loss()