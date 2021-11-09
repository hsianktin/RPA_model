push!(LOAD_PATH,pwd())
using TonksGaswithReactionGaps
simupath = "$(pwd())/data_simu"
function modulate_data(EX,ET) # record the state of data every second.
    Time = [i for i in 0:1:2400]
    Ext = [0.0 for i in Time]
    for i in 1:length(EX)
        t = ceil(Int,ET[i])
        Ext[t+1:end] .= EX[i]
    end
    return [Ext,Time]
end

function modulate_data(Gaps,Time,l1,l2,type = "cumulative")
    T = [i for i in 0:1:2400]
    if type == "cumulative"
        Mod_Gap = Array{Any,1}[]
        for i in 1:l1+l2
            Gap_i = [0 for t in T]
            for j in 1:length(Time)
                t_n = ceil(Int,Time[j])
                Gap_i[t_n+1:end] .= Gaps[i][j]
            end
            push!(Mod_Gap,Gap_i)
        end
        return Mod_Gap
    else
        Mod_Gap = [Array{Int32,1}() for t in T]
        for i in 1:length(Time)
            t_n = ceil(Int,Time[i])
            for j in t_n+1:length(Mod_Gap)
                Mod_Gap[j] = Gaps[i]
            end
        end
        return Mod_Gap
    end
end

Paras = ARGS

l1 = 30
l2 = 20
k_on = parse(Float64,Paras[1])
k_off = parse(Float64,Paras[2])
v_open = parse(Float64,Paras[3])
v_close = parse(Float64,Paras[4])
fold = parse(Float64,Paras[5])
L = parse(Int64,Paras[6])
T1 = parse(Float64,Paras[7])
T2 = parse(Float64,Paras[8])
N = parse(Int64,Paras[9])
exp_label = Paras[10]
if exp_label == "''"
    exp_label = ""
end
simu_label = Paras[11]
gaps_type = Paras[12]
saveid = "$(exp_label)_$(simu_label)"
if exp_label == "wt_15mM_salt_with_diffusion"
    D1 = 3.86
    D2 = 5.22
elseif exp_label == "wt_150mM_salt_with_diffusion"
    D1 = 3.96
    D2 = 4.58
else
    D1 = D2 = 0.0
end

if gaps_type == "none"
    paras = [k_on,k_off,v_open,v_close,fold,L,T1,T2,N]
    paras = [convert(Float64,x) for x in paras]
    function simulation(T1,T2,Time,Length,state_1,state_2,state)
        Simulate(v_open,v_close,0.0,0.0,k_on,k_off,D1,D2,L,l1,l2,T1,Time,Length,state_1,state_2,state) # we assume that the state of binding length 20 is an absolute intermediate state. If ABCD(l=30) wants to detach, or bind, it must go through ABC(l=20) state first. D1 and D2 are assumed similarly.
        Simulate(v_open,v_close,0.0,0.0,k_on*fold,k_off,D1,D2,L,l1,l2,T2,Time,Length,state_1,state_2,state)
    end
    function simu_flow(T1,T2)
        Time = Array{Float64,1}()
        push!(Time,0)
        Length = Array{Float64,1}()
        state_1 = Array{Float64,1}()
        state_2 = Array{Float64,1}()
        # Gaps = Array{Int64,1}[]
        # i = 1
        # while i ≤ l1+l2
        #     push!(Gaps,[L+1-i])
        #     i += 1
        # end
        push!(state_1,0)
        push!(state_2,0)
        push!(Length,0)
        state = Array{Int64,1}()
        for i in 1:L
            push!(state,0)
        end
        simulation(T1,T2,Time,Length,state_1,state_2,state)
        state_1, _ = modulate_data(state_1, Time)
        state_2, _ = modulate_data(state_2, Time)
        # gaps = modulate_data(Gaps,Time,l1,l2)
        storaged_data_line_1 = [paras;[0.0];state_1]'
        storaged_data_line_2 = [paras;[1.0];state_2]'
        f = open("$simupath/rsa_plot_$saveid.csv","a")
            writedlm(f,storaged_data_line_1,",")
            writedlm(f,storaged_data_line_2,",")
        close(f)
    end
    function run_flow(T1,T2,N)
        i = 1
        for  i in  0:N
            if i == 0
            else
                # println()
                simu_flow(T1,T2)
                # if i < N
                    # print("\u1b[1F")
                    # print("\u1b[0K")
                # end
            end
        end
    end
    run_flow(T1,T2,N)
elseif gaps_type == "cumulative"
    paras = [k_on,k_off,v_open,v_close,fold,L,T1,T2,N]
    paras = [convert(Float64,x) for x in paras]
    function simulation(T1,T2,Time,Length,Gaps,state_1,state_2,state)
        Simulate(v_open,v_close,0.0,0.0,k_on,k_off,D1,D2,L,l1,l2,T1,Time,Length,Gaps,state_1,state_2,state,gaps_type) # we assume that the state of binding length 20 is an absolute intermediate state. If ABCD(l=30) wants to detach, or bind, it must go through ABC(l=20) state first. D1 and D2 are assumed similarly.
        Simulate(v_open,v_close,0.0,0.0,k_on*fold,k_off,D1,D2,L,l1,l2,T2,Time,Length,Gaps,state_1,state_2,state,gaps_type)
    end
    function simu_flow(T1,T2)
        Time = Array{Float64,1}()
        push!(Time,0)
        Length = Array{Float64,1}()
        state_1 = Array{Float64,1}()
        state_2 = Array{Float64,1}()
        Gaps = Array{Int64,1}[]
        i = 1
        while i ≤ l1+l2
            push!(Gaps,[L+1-i])
            i += 1
        end
        push!(state_1,0)
        push!(state_2,0)
        push!(Length,0)
        state = Array{Int64,1}()
        for i in 1:L
            push!(state,0)
        end
        simulation(T1,T2,Time,Length,Gaps,state_1,state_2,state)
        state_1, _ = modulate_data(state_1, Time)
        state_2, _ = modulate_data(state_2, Time)
        gaps = modulate_data(Gaps,Time,l1,l2)
        storaged_data_line_1 = [paras;[0.0];state_1]'
        storaged_data_line_2 = [paras;[1.0];state_2]'
        f = open("$simupath/rsa_plot_$saveid.csv","a")
            writedlm(f,storaged_data_line_1,",")
            writedlm(f,storaged_data_line_2,",")
        close(f)
        f = open("$simupath/rsa_cumulative_gaps_$saveid.csv","a")
            for i in 1:(l1+l2)
                writedlm(f,[paras;[i];gaps[i]]',",")
            end
        close(f)
    end
    function run_flow(T1,T2,N)
        i = 1
        for  i in  0:N
            if i == 0
            else
                # println()
                simu_flow(T1,T2)
                # if i < N
                    # print("\u1b[1F")
                    # print("\u1b[0K")
                # end
            end
        end
    end
    run_flow(T1,T2,N)
else
    paras = [k_on,k_off,v_open,v_close,fold,L,T1,T2,N]
    paras = [convert(Float64,x) for x in paras]
    function simulation(T1,T2,Time,Length,Gaps,state_1,state_2,state)
        Simulate(v_open,v_close,0.0,0.0,k_on,k_off,D1,D2,L,l1,l2,T1,Time,Length,Gaps,state_1,state_2,state,gaps_type) # we assume that the state of binding length 20 is an absolute intermediate state. If ABCD(l=30) wants to detach, or bind, it must go through ABC(l=20) state first. D1 and D2 are assumed similarly.
        Simulate(v_open,v_close,0.0,0.0,k_on*fold,k_off,D1,D2,L,l1,l2,T2,Time,Length,Gaps,state_1,state_2,state,gaps_type)
    end
    function simu_flow(T1,T2)
        Time = Array{Float64,1}()
        push!(Time,0)
        Length = Array{Float64,1}()
        state_1 = Array{Float64,1}()
        state_2 = Array{Float64,1}()
        Gaps = Array{Int32,1}[]
        push!(Gaps,[L])
        push!(state_1,0)
        push!(state_2,0)
        push!(Length,0)
        state = Array{Int64,1}()
        for i in 1:L
            push!(state,0)
        end
        simulation(T1,T2,Time,Length,Gaps,state_1,state_2,state)
        state_1, _ = modulate_data(state_1, Time)
        state_2, _ = modulate_data(state_2, Time)
        gaps = modulate_data(Gaps,Time,l1,l2,gaps_type)
        storaged_data_line_1 = [paras;[0.0];state_1]'
        storaged_data_line_2 = [paras;[1.0];state_2]'
        f = open("$simupath/rsa_plot_$saveid.csv","a")
            writedlm(f,storaged_data_line_1,",")
            writedlm(f,storaged_data_line_2,",")
        close(f)
        f = open("$simupath/rsa_exact_gaps_$saveid.csv","a")
            for i in 1:length(gaps)
                writedlm(f,[paras;[i];gaps[i]]',",")
            end
        close(f)
    end
    function run_flow(T1,T2,N)
        i = 1
        for  i in  0:N
            if i == 0
            else
                # println()
                simu_flow(T1,T2)
                # if i < N
                    # print("\u1b[1F")
                    # print("\u1b[0K")
                # end
            end
        end
    end
    run_flow(T1,T2,N)
end
