# Basis for evaluation related functions. Load necessary data and provide a set of functions to evaluate the performance of the parameters.
using Statistics
using CSV
using Plots
using DataFrames
using DelimitedFiles
using Printf
simupath = "$(pwd())/data_simu"
mainpath = pwd()
exppath = "$mainpath/data_exp"
figpath = "$mainpath/figs"
k_ons = Array{Float64,1}()
k_offs = Array{Float64,1}()
v_opens = Array{Float64,1}()
v_closes = Array{Float64,1}()
data_folds = Array{Float64,1}() # used for recording data
Ls = Array{Float64,1}()
T1s = Array{Float64,1}()
T2s = Array{Float64,1}()
Ns = Array{Float64,1}()
state_1s = Array{Float64,1}[]
state_2s = Array{Float64,1}[]
# X_inis = Array{Float64,1}() # experimental values of unprocessed norms
# D_inis = Array{Float64,1}()
intensity_weight = 1

include("intensity_utils.jl")

function dbpath(exp_label,simu_label)
    path = "$simupath/rsa_plot_$(exp_label)_$(simu_label).csv"
    return path
end
function load(dbpath)
    f = open(dbpath,"r")
    n = countlines(f)
    seekstart(f)
    try 
        for i in 1:floor(Int32,n/2)
            # read simu_data.
            temp_line_1 = [parse(Float64,i) for i in split(readline(f),",")]
            temp_line_2 = [parse(Float64,i) for i in split(readline(f),",")]
            # temporal QC 
            # if length(temp_line_1) != 2411
            #     # println(length(temp_line_1))
            #     continue
            # elseif length(temp_line_2) != 2411
            #     # println(length(temp_line_2))
            #     continue
            # end
            k_on,k_off,v_open,v_close,fold,L,T1,T2,N,marker = temp_line_1[1:10]
            if marker == 0.0 # marker is used to ensure sequential order, aka the data is intact
                state_1 = temp_line_1[11:end]
                push!(k_ons,k_on)
                push!(k_offs,k_off)
                push!(v_opens,v_open)
                push!(v_closes,v_close)
                push!(data_folds,fold)
                push!(Ls,L)
                push!(T1s,T1)
                push!(T2s,T2)
                push!(Ns,N)
                push!(state_1s,state_1)
            end
            k_on,k_off,v_open,v_close,fold,L,T1,T2,N,marker = temp_line_2[1:10]
            if marker == 1.0
                state_2 = temp_line_2[11:end]
                push!(state_2s,state_2)
                @printf("loading %s: %.3f\r",dbpath,(i/ceil(Int,n/2)))
            else
                println("input not match")
            end
        end
    catch
    end
    close(f)
    # data_frame = DataFrame(k_on=k_ons,k_off=k_offs,v_open=v_opens,v_close=v_closes,fold=data_folds,L=Ls,state_1=state_1s,state_2=state_2s)
end
function index(k_ons,k_offs,v_opens,v_closes,data_folds,folds)
    index_p = DataFrame(k_on= Float64[], k_off= Float64[], v_open= Float64[], v_close= Float64[])
    # index_p = zeros(0,6)
    for k_on in K_ON
        inds_k_on = findall(x->x==k_on,k_ons)
        for k_off in K_OFF
            inds_k_off = findall(x->x==k_off,k_offs)
            for v_open in V_OPEN
                inds_v_open = findall(x->x==v_open,v_opens)
                for v_close in V_CLOSE
                    inds_v_close = findall(x->x==v_close,v_closes)
                    para_flag = true
                    for fold in folds
                        filter_folds = findall(x->x==fold,data_folds)
                        if length(intersect(inds_k_on,inds_k_off,inds_v_open,inds_v_close,filter_folds)) == 0
                            para_flag = false
                        end
                    end
                    if para_flag
                        push!(index_p,[k_on,k_off,v_open,v_close])
                    end
                end
            end
        end
    end
    # print(length(index_p.k_on))
    return index_p
end
function initialize(exp_label,simu_label)
    load(dbpath(exp_label,simu_label))
    global K_ON = unique(k_ons)
    global K_OFF = unique(k_offs)
    global V_OPEN = unique(v_opens)
    global V_CLOSE = unique(v_closes)
    global FOLDS = unique(data_folds)
    global LS = unique(Ls)
    global T1S = unique(T1s)
    global T2S = unique(T2s)
    global exp_dict = Dict{String,Array{String,1}}()

    if exp_label == "wt_15mM_salt"
        global exp_data_base=[("wt",0),("wt",1),("wt",4),("wt",10),("wt",25),("wt",50)]
        global folds = [0,1,4,10,25]
    elseif exp_label == "wt_150mM_salt"
        global exp_data_base=[("lwt",0),("lwt",1),("lwt",4),("lwt",10),("lwt",25),("lwt",50)]
        global folds = [0,1,4,10,25]
    elseif exp_label == "general" # don't use this now
        global exp_data_base=[("wt",0),("wt",0.1),("wt",0.4),("wt",1.0),("wt",2.5),("lwt",0),("lwt",0.1),("lwt",0.4),("lwt",1.0),("lwt",2.5)]
        global folds = [0,4,10,25]
    else
        global exp_data_base=[("lwt",0),("lwt",0.1),("lwt",0.4),("lwt",1.0),("lwt",2.5)]
        global folds = [0,1,4,10,25]
    end
    function salt_concentration(type) # historically, the wt label refers to wt at 15mM salt; lwt for 150mM salt.
        if type == "wt"
            return "15mM salt"
        elseif type == "lwt"
            return "150mM salt"
        end
    end
    function exp_dict_inject(type::String,concentration::Float64)
        datapath = "$(mainpath)/data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$concentration"]=[datapath,fig_label]
    end
    function exp_dict_inject(type::String,concentration::Int)
        datapath = "$(mainpath)/data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$(convert(Float64,concentration))"]=[datapath,fig_label]
    end
    function exp_dict_inject(tuple::Tuple)
        if length(tuple) == 2
            type,concentration = tuple
            exp_dict_inject(type,concentration)
        else
            println("input not match")
        end
    end
    exp_dict_inject.(exp_data_base)
    # global L = 1000
    global T1 = mean(T1S)
    global T2 = mean(T2S)
    # global l1 = 30
    # global l2 = 20
    global index_p = index(k_ons,k_offs,v_opens,v_closes,data_folds,folds)
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

function trace_statistics(state_1_collection,state_2_collection,α,β) # convert simulation data the mean observed normalized extension
    time_course = [t for t in 1:ceil(T1+T2)]
    function mean_trace(state_collection)
        μ_state = Array{Float64,1}()
        σ_state = Array{Float64,1}()
        time_course = [t for t in 1:ceil(T1+T2)]
        for i in 1:length(time_course)
            X_temp = [x[i] for x in state_collection]
            push!(μ_state,mean(X_temp))
            push!(σ_state,std(X_temp))
        end
        return [μ_state,σ_state]
    end
    μ_state_1,σ_state_1 = mean_trace(state_1_collection)
    μ_state_2,σ_state_2 = mean_trace(state_2_collection)
    function norm(μ_state_1,μ_state_2,α,β)
        μ_1 = mean([μ_state_1[1775],μ_state_1[1790]])
        μ_2 = mean([μ_state_2[1775],μ_state_2[1790]])
        sum = α*μ_1 + β*μ_2 + (1-μ_1-μ_2)
        return sum
    end
    μ_X = [α*μ_state_1[i] + β*μ_state_2[i] + (1-μ_state_1[i]-μ_state_2[i]) for i in eachindex(μ_state_1)]
    μ_X = μ_X./norm(μ_state_1,μ_state_2,α,β)
    σ_X = ((α.*σ_state_1).+(β.*σ_state_2))./norm(μ_state_1,μ_state_2,α,β)
    return (time_course,μ_X,σ_X)
end

function microstates(state_1_collection,state_2_collection) # convert simulation data the mean observed normalized extension
    time_course = [t for t in 1:ceil(T1+T2)]
    function mean_trace(state_collection)
        μ_state = Array{Float64,1}()
        σ_state = Array{Float64,1}()
        time_course = [t for t in 1:ceil(T1+T2)]
        for i in 1:length(time_course)
            X_temp = [x[i] for x in state_collection]
            push!(μ_state,mean(X_temp))
            push!(σ_state,std(X_temp))
        end
        return [μ_state,σ_state]
    end
    μ_state_1,σ_state_1 = mean_trace(state_1_collection)
    μ_state_2,σ_state_2 = mean_trace(state_2_collection)
    return (time_course,μ_state_1,μ_state_2,σ_state_1,σ_state_2)
end

function trace_plot!(state_1_collection,state_2_collection,α,β,label="")
    t_X,μ_X,σ_X = trace_statistics(state_1_collection,state_2_collection,α,β)
    plot!(t_X,μ_X,ribbon=σ_X,fillalpha=.2,lw=5,line=:dash,label=label) 
    yaxis!(:flip)
end

function access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
    DATA=access(k_on,k_off,v_open,v_close,fold)
    state_1_collection = DATA[1]
    state_2_collection = DATA[2]
    return trace_statistics(state_1_collection,state_2_collection,α,β)
end

function access_microstates(k_on,k_off,v_open,v_close,fold)
    DATA=access(k_on,k_off,v_open,v_close,fold)
    state_1_collection = DATA[1]
    state_2_collection = DATA[2]
    return microstates(state_1_collection,state_2_collection)
end

function access_intensity_trace(k_on,k_off,v_open,v_close,fold)
    time,μ_state_1,μ_state_2,σ_state_1,σ_state_2 = access_microstates(k_on,k_off,v_open,v_close,fold)
    intensity = μ_state_1 + μ_state_2
    # normalize intensity when time ∈ (1775,1790)
    relative_intensity = intensity./mean(intensity[ind for ind in eachindex(time) if time[ind] ∈ (1775,1790)])
    return relative_intensity, time
end

function access_trace_statistics(paras::Array{Float64,1})
    if length(paras) == 7
        k_on,k_off,v_open,v_close,fold,α,β = paras
        return access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
    else
        println("length not match")
    end
end

function access_microstates(paras::Array{Float64,1})
    if length(paras) == 5
        k_on,k_off,v_open,v_close,fold = paras
        return access_access_microstates(k_on,k_off,v_open,v_close,fold,α,β)
    else
        println("length not match")
    end
end

function access_trace_statistics(datapath::String) # access experiment data
    readdir(datapath)
    X=Array{Float64,1}[]
    T=Array{Float64,1}[]
    for fname in readdir(datapath)
        # if length(fname) == 12 
            f = open("$datapath/$fname", "r")
            length_file = countlines(f)
            seekstart(f)
            X_temp = [parse(Float64,readline(f)) for i in 1:length_file] # new notation
            T_temp = [5i for i in 1:size(X_temp)[1]]
            push!(X,X_temp[1:481])
            push!(T,T_temp[1:481])
            close(f)
        # end
    end
    ## We should notice that different trace represents DNA of different lengths.
    # X_norm = mean([mean(x[355,358]) for x in X])
    # D_norm = mean([std(x[355:358]) for x in X])
    # push!(X_inis,X_norm)
    # push!(D_inis,D_norm)
    X_backup = X
    X = [x/mean([x[355],x[358]]) for x in X]
    time_course = T[1]
    μ_X = Array{Float64,1}()
    σ_X = Array{Float64,1}()
    for i in 1:size(time_course)[1]
        X_temp = [x[i] for x in X]
        push!(μ_X,mean(X_temp))
        push!(σ_X,std(X_temp))
    end
    # μ_X = μ_X
    # σ_X = σ_X
    # open("$datapath/$fig_label.csv","w") do io
    #     writedlm(io,[time_course μ_X σ_X])
    # end
    X = X_backup
    # plot!(time_course,μ_X,ribbon=σ_X,fillalpha=.1,lw=5,label=fig_label) 
    return (time_course,μ_X,σ_X)
end
function trace_plot!(datapath,fig_label)
    # readdir(datapath)
    # X=Array{Float64,1}[]
    # T=Array{Float64,1}[]
    # for fname in readdir(datapath)
    #     if length(fname) == 12
    #     f = open("$datapath/$fname", "r")
    #     X_temp = [parse(Float64,i) for i in split(readline(f), "," )]
    #     T_temp = [5i for i in 1:size(X_temp)[1]]
    #     push!(X,X_temp[1:481])
    #     push!(T,T_temp[1:481])
    #     close(f)
    #     end
    # end
    # ## We should notice that different trace represents DNA of different lengths.
    # X_norm = mean([mean(x[355,358]) for x in X])
    # D_norm = mean([std(x[355,358]) for x in X])
    # push!(X_inis,X_norm)
    # push!(D_inis,D_norm)
    # X_backup = X
    # X = [x/mean(x[355,358]) for x in X]
    # time_course = T[1]
    # μ_X = Array{Float64,1}()
    # σ_X = Array{Float64,1}()
    # for i in 1:size(time_course)[1]
    #     X_temp = [x[i] for x in X]
    #     push!(μ_X,mean(X_temp))
    #     push!(σ_X,std(X_temp))
    # end
    # μ_X = μ_X
    # σ_X = σ_X
    # open("$datapath/$fig_label.csv","w") do io
    #     writedlm(io,[time_course μ_X σ_X])
    # end
    time_course, μ_X, σ_X = access_trace_statistics(datapath)
    plot!(time_course,μ_X,ribbon=σ_X,fillalpha=.1,lw=5,label=fig_label) 
    return (time_course, μ_X, σ_X, fig_label)
end

function ensemble_plot(k_on,k_off,v_open,v_close,folds,α,β)
    fig=plot(size=(800,1600))
    for fold in folds
        conc = convert(Float64,fold)
        trace_plot!(exp_dict["$conc"][1],exp_dict["$conc"][2])
    end
    for fold in folds
        t_X,μ_X,σ_X=access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
        t_X=[i for i in 1:length(μ_X)]
        temp_df = DataFrame(time=t_X,extension=μ_X,std=σ_X)
        CSV.write("./figs/plot_$(exp_label)_$(simu_label)_$fold.csv",temp_df)
        label = @sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β)
        plot!(t_X,μ_X,line=:dash,lw=5,label = label)
    end
end

function ensemble_plot_int(k_on,k_off,v_open,v_close,simu_folds,exp_folds,α,β,exp_label,simu_label)
    plot_df = DataFrame(T=Float64[], X=Float64[], ΔX = Float64[], label= String[])
    # fig=plot(size=(800,1600))
    for fold in exp_folds
        conc = convert(Float64,fold)
        t_X,μ_X,σ_X,label = trace_plot!(exp_dict["$conc"][1],exp_dict["$conc"][2])
        for i in 1:length(t_X)
            push!(plot_df, [t_X[i],μ_X[i], σ_X[i], label ])
        end
    end
    for fold in simu_folds
        t_X,μ_X,σ_X=access_trace_statistics(k_on,k_off,v_open,v_close,fold,α,β)
        t_X=[i for i in 1:length(μ_X)]
        temp_df = DataFrame(time=t_X,extension=μ_X,std=σ_X)
        CSV.write("./figs/plot_$(exp_label)_$(simu_label)_$fold.csv",temp_df)
        label = "fold = $fold"
        # plot!(t_X,μ_X,line=:dash,lw=5,label = label)
        for i in 1:length(t_X)
            push!(plot_df, [t_X[i],μ_X[i], σ_X[i], label ])
        end
    end
    p=Plotly.plot(
        plot_df,
        x=:T,
        y=:X,
        error_y=:ΔX,
        group=:label,
        title = @sprintf("k_on=%.1E,k_off=%.1E,v_open=%.1E,v_close=%.1E,α=%.2f,β=%.2f",k_on,k_off,v_open,v_close,α,β), 
    )
    CSV.write("./figs/sources/ensemble_plot_df_$(exp_label)_$(simu_label).csv",plot_df)
    open("figs/ensemble_plot_$(exp_label)_$(simu_label).html","w") do io
        PlotlyBase.to_html(io,p.plot)
    end
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
        push!(plot_array,plot(t_X,[μ_1,μ_2],lw=5,ribbon=[σ_1,σ_2],fillalpha=0.2,title="$(label)",legend=false)) # make a plot and add it to the plot_array
    end
    if length(folds)==5
    plot(plot_array...,layout=@layout([a;b;c;d;e]),size=(800,800))
    elseif length(folds)==6
        plot(plot_array...,layout=@layout([a;b;c;d;e;f]),size=(800,800))
    end
    savefig("./figs/microstates_$(exp_label)_$(simu_label).png")
end

function microstates_plot(k_on,k_off,v_open,v_close,folds,exp_label,simu_label)
    plot_df = DataFrame(T=Float64[], mode30nt=Float64[], mode20nt=Float64[], std_mode30nt=Float64[],std_mode20nt=Float64[], label=String[])
    for fold in folds
        t_X,μ_1,μ_2,σ_1,σ_2 = access_microstates(k_on,k_off,v_open,v_close,fold)
        label = @sprintf("fold = %.1f",fold)
        temp_df = DataFrame(time=t_X,ABCD=μ_1,ABC=μ_2,sd_ABCD=σ_1,sd_ABC=σ_2)
        for i in 1:length(t_X)
            push!(plot_df, [t_X[i], μ_1[i], μ_2[i], σ_1[i], σ_2[i], label])
        end
        CSV.write("./figs/microstates_$(exp_label)_$(simu_label)_$fold.csv",temp_df)
        # fig = plot()
        # plot!(t_X,μ_2,lw=5,ribbon=σ_2,fillalpha=0.2,legend=false)
        # push!(plot_array,plot(t_X,[μ_1,μ_2],lw=5,ribbon=[σ_1,σ_2],fillalpha=0.2,title="$(label)",legend=false)) # make a plot and add it to the plot_array
    end
    CSV.write("figs/sources/microstate_$(exp_label)_$(simu_label).csv",plot_df)
    p1 = Plotly.plot(
        plot_df,
        x = :T,
        y = :mode30nt,
        group = :label,
        # error_y = 
    )
    CSV.write("./figs/sources/microstates_plot_df_$(exp_label)_$(simu_label).csv",plot_df)
    open("figs/microstates_plot_1_$(exp_label)_$(simu_label).html","w") do io
        PlotlyBase.to_html(io,p1.plot)
    end
    p2 = Plotly.plot(
        plot_df,
        x = :T,
        y = :mode20nt,
        group = :label,
        # error_y = 
    )
    open("figs/microstates_plot_2_$(exp_label)_$(simu_label).html","w") do io
        PlotlyBase.to_html(io,p2.plot)
    end
    # if length(folds)==5
    # plot(plot_array...,layout=@layout([a;b;c;d;e]),size=(800,800))
    # elseif length(folds)==6
    #     plot(plot_array...,layout=@layout([a;b;c;d;e;f]),size=(800,800))
    # end
    # savefig("./figs/microstates_$(exp_label)_$(simu_label).png")
end

function access_occupancy(k_on,k_off,v_open,v_close,fold)
    t_X,μ_1,μ_2,σ_1,σ_2 = access_microstates(k_on,k_off,v_open,v_close,fold)
    μ = μ_1[end]+μ_2[end]
    σ = √(σ_1[end]^2+σ_2[end]^2)
    return μ,σ
end

function occupancy_plot(k_on,k_off,v_open,v_close,simu_folds,exp_label,simu_label)
    μₒ = Array{Float64,1}()
    σₒ = Array{Float64,1}()
    for fold in simu_folds
        μ_X,σ_X = access_occupancy(k_on,k_off,v_open,v_close,fold)
        push!(μₒ,μ_X)
        push!(σₒ,σ_X)
    end
    plot(simu_folds,μₒ,ribbon=σₒ,fillalpha=.2,line=:dash,lw=5,label = "$(exp_label)")
    savefig("./figs/occupancy_$(exp_label)_$(simu_label).png")
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


function diff(EX,μ_X; type="squared errors")
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
        paras = [k_on,k_off,v_open,v_close,fold,α,β]
        paras = [convert(Float64,x) for x in paras]
        t_X,μ_X,σ_X = access_trace_statistics(paras)
        t_I,μ_I,σ_I = access_intensity_trace(k_on, k_off, v_open, v_close, fold)
        if isnan(μ_X[1]) || isnan(μ_I[1])
            push!(Diffs,100.0)
        else
            entry = exp_dict["$(convert(Float64,fold))"]
            T,X,D = access_trace_statistics(entry[1]) # obtaining the experimental trace
            I, T_I = access_intensity_trace(exp_label, fold)
            EI = [0.0 for i in 1:length(μ_I)]
            EX = [0.0 for i in 1:length(μ_X)] # padding the experimental trace
            # such that the length of the experimental trace is the same as the simulated trace
            # the simulated trace is sampled with dt = 1s, while the experimental trace is sampled with dt = 5s
            # The code works only if the experimental trace is sampled at a larger time interval than the simulated trace
            for i in eachindex(T)
                EX[maximum([ceil(Int,T[i]),1]):end].=X[i]
                EI[maximum([ceil(Int,T_I[i]),1]):end].+=I[i]
            end
            push!(Diffs,Weight[i]*(diff(EX,μ_X;type = "squared errors") + intensity_weight * diff(EI,μ_I;type = "squared errors")))
        end
    end
    return Diffs
end
function evaluate(df,k_on,k_off,v_open,v_close,folds)

    for α in A
        B = [i*α/5 for i in 7:20]
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

function evaluate(df,k_on,k_off,v_open,v_close,folds,α,β)
    difference = sum(diff(k_on,k_off,v_open,v_close,folds,α,β))
    push!(df,[k_on,k_off,v_open,v_close,α,β,difference])
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


# this function is reserved primarily for retraction analysis, which reuses the 'difference' channel in the dataframe.
# function retract(μ_X)
#     μ_X = μ_X[1800:end]
#     max_val,ind=findmax(μ_X)
#     max_val=mean(μ_X[maximum([ind-1,1]):minimum([ind+4,length(μ_X)])])
#     eq_val = mean(μ_X[end-5:end])
#     return max_val-eq_val
# end
# function retraction_watcher(df,k_on_1,k_off,D,L,l,T1,T2,k_off_0,p)
#     retracts = Array{Float32,1}()
#     for i in 1:length(folds)
#         fold = folds[i]
#         k_on_2 = k_on_1*fold
#         paras = [k_on_1,k_on_2,k_off,D,L,l,T1,T2,k_off_0,p]
#         paras = [convert(Float64,x) for x in paras]
#         time_course, μ_X = access_trace_statistics(paras)
#         if isnan(μ_X[1])
#             push!(retracts,0.0)
#         else
#             push!(retracts,retract(μ_X))
#         end
#     end
#     push!(df,[k_on_1,k_off,k_off_0,D,p,0,maximum(retracts)])
# end
# function retract_analyze(df,para)
#     X = sort(unique(df[:,para]))
#     Y = Array{Float64,1}()
#     for x in X
#         temp_df = df[df[:,para].==x,:]
#         push!(Y,maximum(temp_df.diff))
#     end
#     plot(X,Y,line=:dash,lw=5,xlabel=para,ylabel="retraction")
#     if minimum(X) > 0
#         if para !="b"
#             xaxis!(:log)
#         end
#     end
#     savefig("$figpath/retraction_$(para)_$(exp_label)_$(simu_label).svg")
# end