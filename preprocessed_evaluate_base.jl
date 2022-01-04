# after preprocessing the data...
# we can use the following code to train the model

using Statistics
using CSV
using Plots
using DataFrames
using FFTW
using DelimitedFiles
using Printf
norm_debug = true
simupath = "$(pwd())/data_simu"
mainpath = pwd()
exppath = "$mainpath/data_exp"
figpath = "$mainpath/figs"

function initialize(exp_label,simu_label)
    # global exp_dict = Dict{String,Array{String,1}}()
    global folds = [0,1,4,10,25]
    # function salt_concentration(type) # historically, the wt label refers to wt at 15mM salt; lwt for 150mM salt.
    #     if type == "wt"
    #         return "15mM salt"
    #     elseif type == "lwt"
    #         return "150mM salt"
    #     end
    # end
    # function exp_dict_inject(fold)
    #     datapath = "$(mainpath)/data/exp_$(exp_label)_$(fold).csv"
    #     fig_label = "$(exp_label) yRPA 0.01nM to $(0.01*fold)nM"
    #     exp_dict["$concentration"]=[datapath,fig_label]
    # end
    
    # exp_dict_inject.(folds)
    
    global index_p = CSV.read("./data/index_$(exp_label)_$(simu_label).csv",DataFrame)
    global df = DataFrame(k_on = Float64[], k_off = Float64[], v_open = Float64[], v_close = Float64[], α = Float64[], β = Float64[], diff = Float64[])
end
function access(k_on,k_off,v_open,v_close,fold)
    inds_k_on = findall(x->x==k_on,k_ons)
    inds_k_off = findall(x->x==k_off,k_offs)
    inds_v_open = findall(x->x==v_open,v_opens)
    inds_v_close = findall(x->x==v_close,v_closes)
    filter_fold = findall(x->x==fold,data_folds)
    indexes = intersect(inds_k_on,inds_k_off,inds_v_open,inds_v_close,filter_fold)
    function index_to_bit(indexes)
        local L = length(k_ons)
        return [i in indexes for i in 1:1:L]
    end
    bit_index = index_to_bit(indexes)
    return [state_1s[bit_index], state_2s[bit_index]]
end


function access_trace_statistics(exp_label,fold) # access experiment data
    datapath = "./Data/exp_$(exp_label)_$(fold).csv"
    test_df = CSV.read(datapath,DataFrame)
    exp_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_X = Array{Float64,1}[],σ_X = Array{Float64,1}[])
    i = 1
    fold = test_df.folds[i]
    T = [parse(Float64,x) for x in split(test_df.T[1][2:end-1],",")]
    μ_X = [parse(Float64,x) for x in split(test_df.μ_X[1][2:end-1],",")] 
    σ_X = [parse(Float64,x) for x in split(test_df.σ_X[1][2:end-1],",")]
    return T, μ_X, σ_X
end

function access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
    test_df = CSV.read("./data/simu_$(exp_label)_$(simu_label)_$(k_on)_$(k_off)_$(v_open)_$(v_close).csv",DataFrame)
    simu_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_1 = Array{Float64,1}[],σ_1 = Array{Float64,1}[], μ_2 = Array{Float64,1}[],σ_2 = Array{Float64,1}[])
    for i in 1:length(test_df.folds)
        fold′ = test_df.folds[i]
        T = [parse(Float64,x) for x in split(test_df.T[i][2:end-1],",")]
        μ_1 = [parse(Float64,x) for x in split(test_df.μ_1[i][2:end-1],",")] 
        σ_1 = [parse(Float64,x) for x in split(test_df.σ_1[i][2:end-1],",")]
        μ_2 = [parse(Float64,x) for x in split(test_df.μ_2[i][2:end-1],",")]
        σ_2 = [parse(Float64,x) for x in split(test_df.σ_2[i][2:end-1],",")]
        push!(simu_df, [fold′, T, μ_1, σ_1, μ_2, σ_2])   
    end
    j = findfirst(x-> x== fold, simu_df.folds)
    T = simu_df.T[j]
    μ_1 = simu_df.μ_1[j]
    σ_1 = simu_df.σ_1[j]
    μ_2 = simu_df.μ_2[j]
    σ_2 = simu_df.σ_2[j]
    μ_X = α.*μ_1 + β.*μ_2
    σ_X = sqrt.(α.*(σ_1.^2) + β.*(σ_2.^2))
    return T,μ_X,σ_X
end

function access_trace_statistics(paras::Array{Float64,1})
    if length(paras) == 7
        k_on,k_off,v_open,v_close,fold,α,β = paras
        return access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
    else
        println("access_trace_statistics length not match")
    end
end

function Gaussian_derivative(X)
    ℓ = length(X)
    σ = 20
    ∂G(x) = 100/√(2π*σ^2) * exp(- x^2/(2σ^2)) * (-x/(σ^2))
    X∂G = zeros(ℓ)
    for i in 1:ℓ
        X∂G[i] = sum([X[minimum([maximum([j,1]),ℓ])]*∂G(i-j) for j in i-20:1:i+20])
    end
    return X∂G
end


function diff(EX,μ_X,type="derivatives")
    if type == "squared errors"
        difference = 0.0
        for j in 1200:2400
            if (j  < 2000.0)
                difference += 1*(EX[j]-μ_X[j])^2
            else
                difference += 2*(EX[j]-μ_X[j])^2
            end
        end
        return difference
    elseif type == "derivatives"
        derivative = 0.0
        for j in 1200:2400
            if (j  < 2000.0)
                derivative += 1*(EX[j]-μ_X[j])^2
            else
                derivative += 2*(EX[j]-μ_X[j])^2
            end
        end
        EX′ = Gaussian_derivative(EX[1800:2400])
        μ_X′ = Gaussian_derivative(μ_X[1800:2400])        
        if norm_debug == true
            println("0-norm: $(derivative), 1-norm: $(sum((EX′.-μ_X′).^2)*5)")
        end
        return derivative + sum((EX′.-μ_X′).^2)*0
    elseif type == "maximum"
        difference = maximum(abs.(EX[1800:2400].-μ_X[1800:2400]))
        return difference
    elseif type == "mixed"
        difference = diff(EX, μ_X, "squared errors")/10 + diff(EX,μ_X)
    end
end

function diff(k_on,k_off,v_open,v_close,folds,α,β)
    Diffs = Array{Float32,1}()
    Weight = [1,1,1,1,1,1,1]
    for i in 1:length(folds)
        fold = folds[i]
        if norm_debug== true
            println("fold=$fold")
        end
        paras = [k_on,k_off,v_open,v_close,fold,α,β]
        paras = [convert(Float64,x) for x in paras]
        t_X,μ_X,σ_X = access_trace_statistics(paras)
        if isnan(μ_X[1])
            push!(Diffs,1.0)
        else
            T,X,D = access_trace_statistics(exp_label,fold)
            EX = [0.0 for i in 1:length(μ_X)]
            for i in 1:length(T)
                EX[maximum([ceil(Int,T[i]),1]):end].=X[i]
            end
            push!(Diffs,Weight[i]*diff(EX,μ_X,"derivatives"))
        end
    end
    return Diffs
end
function evaluate(df,k_on,k_off,v_open,v_close,folds)

    for α in A
        B = [i*α/5 for i in 5:20]
        for β in B
        # β = 2*α
            difference = sum(diff(k_on,k_off,v_open,v_close,folds,α,β))
            push!(df,[k_on,k_off,v_open,v_close,α,β,difference])
            # println([k_on,k_off,v_open,v_close,α,β,difference])
        end
        # only enabled in lab 
        # ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β);
        # savefig("D:\\rpa\\figs\\plot_$(k_on)_$(k_off)_$(v_open)_$(v_close)_$(α)_$(β).png")
    end 
end

function analyze(df,para)
    X = sort(unique(df[:,para]))
    Y = Array{Float64,1}()
    for x in X
        temp_df = df[df[:,para].==x,:]
        push!(Y,minimum(temp_df.diff))
    end
    plot(X,Y,line=:dash,lw=5,xlabel=para,ylabel="min_diff")
    if minimum(X) > 0
        if para !="b"
            xaxis!(:log)
        end
    end
    temp_df = DataFrame(para=X, error=Y)
    CSV.write("$figpath/sources/landscape_$(para)_$(exp_label)_$(simu_label).csv",temp_df)
    savefig("$figpath/landscape/landscape_$(para)_$(exp_label)_$(simu_label).svg")
end


function trace_plot!(exp_label,fold)
    time_course, μ_X, σ_X = access_trace_statistics(exp_label,fold)
    plot!(time_course,μ_X,ribbon=σ_X,fillalpha=.1,lw=5) 
end

function ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β)
    fig=plot(size=(800,600))
    for fold in folds
        conc = convert(Float64,fold)
        trace_plot!(exp_label,fold)
    end
    for fold in folds
        t_X,μ_X,σ_X=access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
        t_X=[i for i in 1:length(μ_X)]
        temp_df = DataFrame(time=t_X,extension=μ_X,std=σ_X)
        CSV.write("./figs/sources/plot_$(exp_label)_$(simu_label)_$fold.csv",temp_df)
        label = @sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β)
        plot!(t_X,μ_X,line=:dash,lw=5,label = label)
    end
end

function access_microstates(k_on,k_off,v_open,v_close,fold)
    test_df = CSV.read("./data/simu_$(exp_label)_$(simu_label)_$(k_on)_$(k_off)_$(v_open)_$(v_close).csv",DataFrame)
    simu_df = DataFrame(folds=Int[],T=Array{Float64,1}[],μ_1 = Array{Float64,1}[],σ_1 = Array{Float64,1}[], μ_2 = Array{Float64,1}[],σ_2 = Array{Float64,1}[])
    for i in 1:length(test_df.folds)
        T = [parse(Float64,x) for x in split(test_df.T[1][2:end-1],",")]
        μ_1 = [parse(Float64,x) for x in split(test_df.μ_1[1][2:end-1],",")] 
        σ_1 = [parse(Float64,x) for x in split(test_df.σ_1[1][2:end-1],",")]
        μ_2 = [parse(Float64,x) for x in split(test_df.μ_2[1][2:end-1],",")]
        σ_2 = [parse(Float64,x) for x in split(test_df.σ_2[1][2:end-1],",")]
        push!(simu_df, [fold, T, μ_1, σ_1, μ_2, σ_2])   
    end
    j = findfirst(x-> x== fold, simu_df.folds)
    T = simu_df.T[j]
    μ_1 = simu_df.μ_1[j]
    σ_1 = simu_df.σ_1[j]
    μ_2 = simu_df.μ_2[j]
    σ_2 = simu_df.σ_2[j]
    return T,μ_1,μ_2,σ_1,σ_2 
end

function microstates_plot(k_on,k_off,v_open,v_close,folds)
    plot_array = Any[] # can type this more strictly
    for fold in folds
        t_X,μ_1,μ_2,σ_1,σ_2 = access_microstates(k_on,k_off,v_open,v_close,fold)
        label = @sprintf("fold = %.1f",fold)
        temp_df = DataFrame(time=t_X,ABCD=μ_1,ABC=μ_2,sd_ABCD=σ_1,sd_ABC=σ_2)
        CSV.write("./figs/microstates_$(exp_label)_$(simu_label)_$fold.csv",temp_df)
        # fig = plot()
        # plot!(t_X,μ_2,lw=5,ribbon=σ_2,fillalpha=0.2,legend=false)
        push!(plot_array,plot(t_X,[μ_1,μ_2],lw=5,ribbon=[σ_1,σ_2],fillalpha=0.2,label=["30nt" "20nt"])) # make a plot and add it to the plot_array
    end
    if length(folds)==5
    plot(plot_array...,layout=@layout([a;b;c;d;e]),size=(800,800))
    elseif length(folds)==6
        plot(plot_array...,layout=@layout([a;b;c;d;e;f]),size=(960,800))
    elseif length(folds)==7
        plot(plot_array...,layout=@layout([a;b;c;d;e;f;g]),size=(1120,800))
    end
    savefig("./figs/microstates_$(exp_label)_$(simu_label).png")
end