# it seems that because frequent exchange of data between CPU and GPU,
# the training process is very slow.
# prefer to use CPU for training and prediction
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
device = gpu
CUDA.allowscalar(false)
# set allow scalar indexing to be false 
# to avoid the error of scalar indexing on GPU


# load data 
df = CSV.read("data_simu/rsa_plot_landscape_mean_L_5000.csv", DataFrame)
# # for the sake of simplicity, we first consider 1D problem
df = df[[row for row in 1:nrow(df) if df.Id[row] == 1], :]
para_df = df[:, 1:5]
# # we also select the data from t = 0 to t = 1800, 
# # to avoid time switching effect temporarily
data_df = df[:, 11:end-120]
# design the neural ODE model
# let x(t) denote the data, p denote the parameters 
# and f(x(t), p) denote the dynamics of x(t)
# dx/dt = f(x(t), p)
# we consider some known information,
#  f(x,p) = p[1] * (1 -x ) - (p[2] + p[4])x + NN(θ, p, x)
# for each row of data_df, we have 
# ts = [0:5:1800] |> collect 
# xs = data_df[row, :] |> collect

# Define the neural network (NN), which takes the current state and parameters as input
NN = Chain(Dense(6, 64, tanh), Dense(64, 1)) |> gpu
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
    xs = [u[1] for u in us]
    return xs
end

# train data = [(get_u0(row), get_xs(row)) for row in 1:nrow(data_df)]
function get_u0(row)
    x0 = data_df[row, 1]
    p0 = para_df[row, :] |> collect
    u0 = vcat(x0, p0)
    return u0
end

function get_xs(row)
    xs = data_df[row, :] |> collect
    return xs
end

train_data = [(get_u0(row), get_xs(row)) for row in 1:nrow(data_df)]  
# convert data to Float32
train_data = [(Float32.(u0), Float32.(xs)) for (u0, xs) in train_data]
# Define the loss function as the sum of squared errors between the true and predicted values
loss(u0, xs) = sum(abs, predict(u0) .- xs) / length(xs)

# parallel computing of du╱dt()
function dU╱dt(U, θ, t)
    # get rows of U as X, Kon, Koff, Kopen, Kclose, Fold
    X, Kon, Koff, Kopen, Kclose, Fold = eachrow(U) 
    known = Kon .* (1 .- X) - (Koff .+ Kclose) .* X
    # dX╱dt = Kon .* (1 .- X) .- (Koff .+ Kclose) .* X .+ ( (re(θ)(
    #     hcat(X, log.(Kon), log.(Koff), log.(Kopen), log.(Kclose), Fold)|> transpose  ) ) |> transpose) .* X .* α
    dX╱dt = known + transpose(re(θ)(
            hcat(X, log.(Kon), log.(Koff), log.(Kopen), log.(Kclose), Fold)|> transpose  ) )  .* X .* α
    dP╱dt = zeros(Float32, (5, length(X))) |> gpu
    ∂ₜU = vcat(dX╱dt |> transpose, dP╱dt)
    return (∂ₜU |> Array |> gpu)
end

U0 = hcat([u0 for (u0, xs) in train_data]...) |> gpu # 6 × N_sample Matrix
Xs = hcat([xs for (u0, xs) in train_data]...) |> gpu # 361 × N_sample Matrix

# U0 = rand(Float32, (6, 2)) |> gpu
# Xs = rand(Float32, (361, 2)) |> gpu


# dU╱dt(U, θ, 0.0) |> typeof # LinearAlgebra.Transpose{Float64, gpu{Float64, 2, CUDA.Mem.DeviceBuffer}}
# typeof(U)
# U + dU╱dt(U, θ, 0.0) .* 0.01

function batch_predict(U, θ)
    Prob = ODEProblem(dU╱dt, U, (0.0f0, 1800.0f0), θ)
    Sol = solve(Prob, Tsit5(), saveat=0:5:1800)
    # Sol.u # 361 Array, convert to (361, 100) Matrix 
    SolU = Sol.u
    return hcat([u[1,:] for u in SolU]...) |> transpose |> gpu
end

batch_predict(U0, θ) |> typeof # 361 × N_sample Matrix

batch_loss(U0, θ, Xs) = sum(abs, batch_predict(U0, θ) .- Xs) / length(Xs)
α = 0.01f0 |> gpu
@time batch_loss(U0, θ, Xs)

# test gradient of loss function
@time ∇θ = gradient(Flux.params(θ)) do
    batch_loss(U0, θ, Xs)
end



# parallel computing of loss() 
function loss()
    l = 0.0
    Threads.@threads for (u0, xs) in train_data
        l += loss(u0, xs)
    end
    l = l / length(train_data)
    return l
end

α = 0 |> device
@info "mean-field ODE loss" loss()

α = 0.001 |> device
@info "untrained neural ODE loss with α = $α" loss()
@time loss()
@time     ∇θ = gradient(Flux.params(θ)) do
    loss()
end
# train 
opt = ADAMW(0.01)

# rearrange train_data by random permutation
# load packages for shuffle
using Random
train_data = train_data[shuffle(1:nrow(data_df)), :]

using ProgressBars
n_epochs = 1000
p = ProgressBar(1:n_epochs)
for iter in p
    # iter = 1
    ∇θ = gradient(Flux.params(θ)) do
        loss()
    end
    Flux.update!(opt, Flux.params(θ), ∇θ)
        set_multiline_postfix(p, "loss = $(loss() |> cpu)")
end
    
@info "trained neural ODE loss with α = $α" loss()

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