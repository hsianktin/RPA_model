# simulation similar to that in FRET experiments

using Random
using ProgressMeter
using CSV
using DataFrames
# using Plots
using QuadGK
# gr()
bin_width = 0.1
## initialize state space
function simu_data_acquire(h_on,h_off,D,threshold)
    free = 30
    l = 30
    Hairpin = 7
    state = zeros(free+Hairpin+3+Hairpin)
    # h_on = 0.1
    # h_off = 0.1
    state[free+1:free+Hairpin] .= -1
    state[1:l] .= 1
    # D = 3
    ## Define function for hairpin dynamics
    function v_h_on(state)
        tmp = 0
        for i in free+1:free+Hairpin
            if state[i] == 0
                tmp = i
            end
        end
        if tmp != 0
            return h_on
        else
            return 0
        end
    end

    function H_ON(state)
        # find the left most 
        tmp = 0
        for i in free+1:free+Hairpin
            if state[i] == 0
                tmp = i
            end
        end
        state[tmp] = -1
    end

    function h_position(state)
        tmp = free+Hairpin
        for i in free+1:free+Hairpin
            if state[i] == -1
                tmp = i
                break
            end
        end
        return tmp
    end

    function v_h_off(state)
        tmp = 0
        for i in free+1:free+Hairpin
            if state[i] == -1
                tmp = i
                break
            end
        end
        if tmp == 0
            return 0 
        else
            return h_off
        end
    end

    function H_OFF(state)
        tmp = 0
        for i in free+1:free+Hairpin
            if state[i] == -1
                tmp = i
                break
            end
        end
        state[tmp] = 0
    end

    function I(i,n) # 从第i个nt开始向右结合的长度为n的蛋白所占据的位置坐标
        [j for j in i:i+n-1]
    end
    ## Cite from LoadOnChains
    function find_protein_position(n,state) #记录蛋白质的位置，确认他们是否相邻，不储存，直接根据state计算
        temp = 0
        temp_dict = Array{Int64,1}()
        i = 1
        while i + n - 1 < size(state)[1] + 1 # 记录的实际上是首位置
            if sum(state[I(i,n)]) == n
                temp += 1
                push!(temp_dict,i)
                i += n # 跳过找到的蛋白质的剩下部分
            else
                i += 1
            end
        end
        return temp_dict 
    end

    function v_diffuse(state)
        L = length(state)
        protein_position = find_protein_position(l,state)
        temp = 0
        for i in 1:size(protein_position)[1] # 相邻的蛋白质不可向左向右diffuse，边缘的不可diffuse
            if protein_position[i] != L-l+1
                if state[protein_position[i]+l] == 0
                    temp+=1
                end
            end
            if protein_position[i] !=1
                if state[protein_position[i]-1] == 0
                    temp+=1
                end
            end
        end
        return temp*D
    end

    function Diffuse(state)
        L = length(state)
        protein_position_1 = find_protein_position(l,state)
        left_move = Array{Int64,1}() # 储存可以向左的蛋白序号（不是位置，是在protein_position数组中的顺序）
        right_move = Array{Int64,1}() # 储存可以向右的序号
        for i in 1:size(protein_position_1)[1]
            if protein_position_1[i] != L-l+1
                if state[protein_position_1[i]+l] == 0
                    push!(right_move,i)
                end
            end
            if protein_position_1[i] !=1
                if state[protein_position_1[i]-1] == 0
                    push!(left_move,i)
                end
            end
        end
        Length = size(left_move)[1] + size(right_move)[1] # 总长度
        j = rand(1:Length)
        if j > size(left_move)[1] # 向右
            j = j - size(left_move)[1] # 获得指定编号
            state[I(protein_position_1[right_move[j]],l)] = zeros(l)
            state[I(protein_position_1[right_move[j]]+1,l)] = ones(l)
        else
            state[I(protein_position_1[left_move[j]],l)] = zeros(l)
            state[I(protein_position_1[left_move[j]]-1,l)] = ones(l)
        end
    end


    function random_select(a,b,c)
        aa = a/(a+b+c)
        bb = (a+b)/(a+b+c)
        cc = (a+b+c)/(a+b+c)
        k = rand()
        if k<aa
            Diffuse(state)
            return randexp(1)[1]/(cc)
        elseif k<bb
            H_ON(state)
            return randexp(1)[1]/(cc)
        else
            H_OFF(state)
            return randexp(1)[1]/cc
        end
    end

    function one_step_firing(state)
        v = v_h_off(state)+v_h_on(state)+v_diffuse(state)
        δt = random_select(v_diffuse(state),v_h_on(state),v_h_off(state))/v
        return δt
    end 

    function sample(state)
        count = 0
        t = 0
        time = Array{Float64,1}()
        push!(time,t)
        readout = Array{Float64,1}()
        #debug
        hairpin_open = Array{Int8,1}()
        push!(readout,maximum([find_protein_position(l,state)[1]+(l-1-free)]))
        push!(hairpin_open,h_position(state))
        is_high_fret=true
        progress_print = 0
        while count<10000
            if is_high_fret
                if find_protein_position(l,state)[1]+(l-1-free) ≥ Hairpin
                    is_high_fret = false 
                end
            else
                if find_protein_position(l,state)[1]+(l-1-free) < Hairpin
                    is_high_fret = true
                    count+=1
                    progress_print=count/10000
                    if floor(Int,progress_print*100)!=floor(Int,(progress_print-1)*100)
                        print("sampling progress: $(floor(Int,progress_print*100)) %\r")
                    end
                end
            end
            t+=one_step_firing(state)
            push!(time,t)
            push!(readout,maximum([find_protein_position(l,state)[1]+(l-1-free)]))
            push!(hairpin_open,h_position(state))
        end
        return [time,readout,hairpin_open]
    end

    time,readout,hairpin_open = sample(state)
    # histogram(readout)

    # statistics of waiting time
    function dwell_time(time,readout,threshold)
        filtered_readout = readout.≥ threshold + 1
        # println(sum(filtered_readout)/length(filtered_readout))
        waiting_time_low_fret = Array{Float64,1}()
        waiting_time_high_fret = Array{Float64,1}()
        time_ini = .0
        state_ini = false
        for i in 2:length(time)
            if filtered_readout[i] != state_ini
                if state_ini
                    push!(waiting_time_low_fret,time[i]-time_ini)
                else 
                    push!(waiting_time_high_fret,time[i]-time_ini)
                end
                state_ini = filtered_readout[i]
                time_ini = time[i]
            end
        end
        return [waiting_time_low_fret,waiting_time_high_fret]
    end
    waiting_time_low_fret,waiting_time_high_fret = dwell_time(time,readout,threshold);
    return waiting_time_low_fret
end
default(legend=true)

function exp_data_acquire(type,salt_concentration::Int)
    path = "./data_exp/fret/exp_data_$(type)_$(salt_concentration).txt"
    f = open(path,"r")
    n = countlines(f)
    seekstart(f)
    low_fret_times = Array{Float64,1}()
    for i in 1:n
        data = parse(Float64,readline(f))
        push!(low_fret_times,data*32/1000)
    end
    return low_fret_times
end

function fixed_hist!(waiting_time_low_fret,label)
    waiting_time_low_fret=waiting_time_low_fret./bin_width
    right_limits = [i for i in 1:1:ceil(Int,maximum(waiting_time_low_fret))]
    values = [0 for i in 1:1:ceil(Int,maximum(waiting_time_low_fret))]
    for value in waiting_time_low_fret
        posi = 1
        while right_limits[posi] < value
            posi += 1
        end
        values[posi]+=1
    end
    bar!((right_limits.-0.5).*bin_width,values./(sum(values)*bin_width),label=label,alpha=0.2)
end

function fixed_hist!(waiting_time_low_fret,D,h_on,h_off,threshold)
    waiting_time_low_fret=waiting_time_low_fret./bin_width
    right_limits = [i for i in 1:1:ceil(Int,maximum(waiting_time_low_fret))]
    values = [0 for i in 1:1:ceil(Int,maximum(waiting_time_low_fret))]
    for value in waiting_time_low_fret
        posi = 1
        while right_limits[posi] < value
            posi += 1
        end
        values[posi]+=1
    end
    bar!((right_limits.-0.5).*bin_width,values./(sum(values)*bin_width),label="simulation D=$D, h_on=$h_on, h_off=$h_off, threshold=$threshold",alpha=0.2)
end

function histogram_plot(waiting_time_low_fret,D,h_on,h_off,threshold)
    t1 = 0.069
    v1 = 1/t1
    t2 = 0.39
    v2 = 1/t2
    f1 = 0.45
    f2 = 0.55
    xs = [i for i in 0:0.002:3]
    ys = [f2*exp(-v2*i)+f1*exp(-v1*i) for i in 0:0.002:3]
    ys = ys./(sum(ys)*0.002)
    xlabel!("dwell time of low fret signals/s")
    ylabel!("mean probability density")
    fixed_hist!(waiting_time_low_fret[waiting_time_low_fret.<3],D,h_on,h_off,threshold)
    # plot!(xs,ys,label="experimental data fits",line=:dash,lw=5)
    xlims!(0,3)
end

function empirical_distribution_on_grids(dwell_time,bin_width,span)
    dwell_time=dwell_time[dwell_time.<span]
    dwell_time=dwell_time./bin_width
    right_limits = [i for i in 1:1:ceil(Int,(span/bin_width))]
    values = [0 for i in 1:1:ceil(Int,(span/bin_width))]
    for value in dwell_time
        posi = 1
        while right_limits[posi] < value
            posi += 1
        end
        values[posi]+=1
    end
    values=values./(sum(values.*bin_width))
    return values
end

function relative_entropy_on_grids(dwell_time_theory,dwell_time_exp,span)
    p = empirical_distribution_on_grids(dwell_time_theory,bin_width,span)
    q = empirical_distribution_on_grids(dwell_time_exp,bin_width,span)
    if sum(p==0) >0 || sum(q==0) >0
        println("error: lack of samples")
        return  Base.Inf64
    else
        S = 0
        for (pᵢ,qᵢ) in zip(p,q)
            S += -qᵢ*log(pᵢ/qᵢ)*bin_width
        end
        return S
    end
end
exp_data_base=[("wt",15),("wt",150),("dwh",15),("dwh",150)]
function comparing_results(type,salt_concentration)
    exp_data = exp_data_acquire(type,salt_concentration)
    h_on = 0.1
    h_off = 0.1
    threshold = 7
    Ts = Array{String,1}()
    Cs = Array{Float64,1}()
    Ss = Array{Float64,1}()
    Ds = [i/10 for i in 1:100]
    D_record = Array{Float64,1}()
    for D in Ds
        dwell_time_low_fret = simu_data_acquire(h_on,h_off,D,threshold)
        push!(Ts,type)
        push!(Cs,salt_concentration)
        push!(Ss,relative_entropy_on_grids(dwell_time_low_fret[dwell_time_low_fret.<maximum(exp_data)],exp_data,3))
        push!(D_record,D)
    end
    bar(D_record,Ss,label="Simulation h_on=$h_on, h_off=$h_off, threshold=$threshold")
    df = DataFrame(Type=Ts,Salt_Concentration=Cs,Diffusion=D_record,Entropy=Ss)
    CSV.write("./data_simu/fret/cross_entropy_$(type)_$(salt_concentration)_$(threshold)_$(bin_width).csv",df)
    xlabel!("D/s^(-1)")
    ylabel!("Cross Entropy against Experiments") 
    savefig("$(pwd())/figs/fret/plot_fret_KLDivergence_D_$(type)_$(salt_concentration)_$(threshold)_$(bin_width).svg")
end

# obtain simulation data
# @showprogress 1 "running" for i in 1:length(exp_data_base)
#     comparing_results(exp_data_base[i][1],exp_data_base[i][2])
# end

# create distribution sample
fitting_parameters = [1.2,2.3,1.1,4.0]
for i in 1:length(exp_data_base)
    exp_data = exp_data_acquire(exp_data_base[i][1],exp_data_base[i][2])
    h_on = 0.1
    h_off = 0.1
    threshold = 7
    dwell_time_low_fret = simu_data_acquire(h_on,h_off,fitting_parameters[i],threshold)
    simu_data = dwell_time_low_fret[dwell_time_low_fret.<maximum(exp_data)]
    plot()
    fixed_hist!(simu_data,fitting_parameters[i],h_on,h_off,threshold)
    fixed_hist!(exp_data,"experiment result of $(exp_data_base[i][2])mM KCl and $(exp_data_base[i][1]) RPA") 
    savefig("$(pwd())/figs/hist_$(exp_data_base[i][1])_$(exp_data_base[i][2]).svg")
end