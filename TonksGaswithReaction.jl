module TonksGaswithReaction
using Random
using DelimitedFiles
export DelimitedFiles, readdlm, writedlm

# state 是长度为L的DNA，0表示未结合蛋白质，1表示已结合
# by default l1 = 30, l2 = 20

function I(i,n) # 从第i个nt开始向右结合的长度为n的蛋白所占据的位置坐标
    [j for j in i:i+n-1]
end

function gaps_count(state,x) # count the number of gaps of size x
    L = length(state)
    count_1 = 0
    for i = 1:L-x+1
        if sum([abs(y) for y in state[i:i+x-1]]) == 0
            count_1 += 1
        end
    end
    return count_1
end

function C(n, state,L) # 记录的是长度为n的所有可结合位点中，可供结合的位点数目
    temp  = 0
    for i in 1:L-n+1
        if sum(state[I(i,n)]) == 0
            temp += 1
        end
    end
    return temp
end

function v_on(state,k_on,l,L) # 计算v_on
    k_on * C(l,state,L)
end

function Lift_1(n, state, L) # 随机结合一个蛋白质
    j = rand(1:C(n,state,L)) # 选取一个位点
    I_D = [k for k in 1:n] # 初始化长度为n的选框，由于julia的语言设定所必须的
    temp  = 0
    for i in 1:L-n+1 # 循环地找第j个位置，因为并没有直接记录gap的位置信息。
        if sum(state[I(i,n)]) == 0 # 找到了空位计数加一
            temp += 1
        end
        if temp == j #找到对应位置，跳出循环，选框赋值
            I_D = I(i,n)
            break
        end
    end
    state[I_D] = ones(n) #在选框的位置都设置为1
end

function Lift_2(n, state, L) # 随机结合一个蛋白质
    j = rand(1:C(n,state,L)) # 选取一个位点
    I_D = [k for k in 1:n] # 初始化长度为n的选框，由于julia的语言设定所必须的
    temp  = 0
    for i in 1:L-n+1 # 循环地找第j个位置，因为并没有直接记录gap的位置信息。
        if sum(state[I(i,n)]) == 0 # 找到了空位计数加一
            temp += 1
        end
        if temp == j #找到对应位置，跳出循环，选框赋值
            I_D = I(i,n)
            break
        end
    end
    state[I_D] = -ones(n) #在选框的位置都设置为-1
end

function find_protein_position_1(n,state) #记录蛋白质的位置，确认他们是否相邻，不储存，直接根据state计算
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

function find_protein_position_2(n,state) #记录蛋白质的位置，确认他们是否相邻，不储存，直接根据state计算
    temp = 0
    temp_dict = Array{Int64,1}()
    i = 1
    while i + n - 1 < size(state)[1] + 1 # 记录的实际上是首位置
        if sum(state[I(i,n)]) == -n
            temp += 1
            push!(temp_dict,i)
            i += n # 跳过找到的蛋白质的剩下部分
        else
            i += 1
        end
        
    end
    return temp_dict 
end

function v_detach(protein_position,k_off)
    return size(protein_position)[1] * k_off
end

function Detach(protein_position,state,l)
    Total = size(protein_position)[1]
    rd = rand(1:Total)
    state[I(protein_position[rd],l)] = zeros(l)
end

function transition_from_1_to_2(protein_position_1,l1,l2,state)
    temp_selector=rand(1:length(protein_position_1))
    state[I(protein_position_1[temp_selector],l2)].=-1
    state[I(protein_position_1[temp_selector]+l2,l1-l2)] .=0
end

function v_transition_from_1_to_2(protein_position_1,l1,l2,state,v_open)
    return v_open*length(protein_position_1)
end

function transition_from_2_to_1(protein_position_2,l1,l2,state)
    lut_available = Array{Int,1}()
    for i in protein_position_2
        if i+l1-1 <= length(state)
            if state[I(i+l2,l1-l2)] == zeros(l1-l2)
                push!(lut_available,i)
            end
        end
    end
    temp_selector=rand(1:length(lut_available))
    state[I(lut_available[temp_selector],l1)].=1
end

function v_transition_from_2_to_1(protein_position_2,l1,l2,state,v_close)
    lut_available = Array{Int,1}()
    for i in protein_position_2
        if i+l1-1 <= length(state)
            if state[I(i+l2,l1-l2)] == zeros(l1-l2)
                push!(lut_available,i)
            end
        end
    end
    return v_close*length(lut_available)
end

function v_diffuse(protein_position,L,l,D,state)
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


function Diffuse_1(protein_position_1,state,L,l)
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

function Diffuse_2(protein_position_2,state,L,l)
    left_move = Array{Int64,1}() # 储存可以向左的蛋白序号（不是位置，是在protein_position数组中的顺序）
    right_move = Array{Int64,1}() # 储存可以向右的序号
    for i in 1:size(protein_position_2)[1]
        if protein_position_2[i] != L-l+1
            if state[protein_position_2[i]+l] == 0
                push!(right_move,i)
            end
        end
        if protein_position_2[i] !=1
            if state[protein_position_2[i]-1] == 0
                push!(left_move,i)
            end
        end
    end
    Length = size(left_move)[1] + size(right_move)[1] # 总长度
    j = rand(1:Length)
    if j > size(left_move)[1] # 向右
        j = j - size(left_move)[1] # 获得指定编号
        state[I(protein_position_2[right_move[j]],l)] = zeros(l)
        state[I(protein_position_2[right_move[j]]+1,l)] = -ones(l)
    else
        state[I(protein_position_2[left_move[j]],l)] = zeros(l)
        state[I(protein_position_2[left_move[j]]-1,l)] = -ones(l)
    end
end

function Diffuse_Kick_Rate(k_off_0,D,l,protein_position)
    # record pairs of adjacent protein
    Pairs = Array{Int64,1}[]
    i = 1
    while i < size(protein_position)[1]
        if protein_position[i+1] - protein_position[i] == l
            push!(Pairs, [i,i+1])
        end
        i+=1
    end
    # facilitated rate
    v_dk = k_off_0*2/3
    return v_dk*size(Pairs)[1]
end

function Diffuse_Kick(l,protein_position,state)
    # record pairs of adjacent protein
    Pairs = Array{Int64,1}[]
    i = 1
    while i < size(protein_position)[1]
        if protein_position[i+1] - protein_position[i] == l
            push!(Pairs, [i,i+1])
        end
        i+=1
    end
    # facilitated rate
    d1 = rand(1:size(Pairs)[1])
    state[I(protein_position[Pairs[d1][1]],l)] = zeros(l)    
    state[I(protein_position[Pairs[d1][2]],l)] = zeros(l)
    state[I(floor(Int,(protein_position[Pairs[d1][1]]+protein_position[Pairs[d1][2]])/2),l)] = ones(l)
end

function char_select(i)
    if i == 1
        return "="
    elseif i == 0
        return "-"
    end
end

function Simulate(v_open,v_close,k_on_1,k_off_1,k_on_2,k_off_2,D1,D2,L,l1,l2,ΔT,Time,Length,Gaps,state_1,state_2,state)
    t = Time[end]
    Δt = 0
    rough_progress_0 = 0
    while Δt < ΔT
        protein_position_1 = find_protein_position_1(l1,state)
        protein_position_2 = find_protein_position_2(l2,state)
        V_on_1 = v_on(state,k_on_1,l1,L)
        V_on_2 = v_on(state,k_on_2,l2,L)
        V_sc_1 = v_transition_from_1_to_2(protein_position_1,l1,l2,state,v_open)
        V_sc_2 = v_transition_from_2_to_1(protein_position_2,l1,l2,state,v_close)
        V_off_1 = v_detach(protein_position_1, k_off_1)
        V_off_2 = v_detach(protein_position_2, k_off_2)
        V_d_1 = v_diffuse(protein_position_1,L,l1,D1,state)
        V_d_2 = v_diffuse(protein_position_2,L,l2,D2,state)
        V_sc = V_sc_1 + V_sc_2
        V_on = V_on_1 + V_on_2
        V_d = V_d_1 + V_d_2
        V_off = V_off_1 + V_off_2
        v = V_on + V_d + V_off + V_sc
        if v == 0
            break
        end
        on_off_over_diffuse = V_d / (V_on + V_d + V_off + V_sc)
        on_over_off = (V_d+V_off) / (V_on + V_d + V_off + V_sc)
        sc_over_on = (V_on + V_d + V_off)/ (V_on + V_d + V_off + V_sc)
        d1 = rand()
        δt = randexp(1)[1]/(V_on + V_d + V_off + V_sc)
        if d1 > sc_over_on
            d2 = rand()
            if d2 > V_sc_1/(V_sc)
                transition_from_2_to_1(protein_position_2,l1,l2,state)
            else
                transition_from_1_to_2(protein_position_1,l1,l2,state)
            end
        elseif d1 > on_over_off
            if V_on_1 == 0
                Lift_2(l2,state,L)
            elseif V_on_2 == 0
                Lift_1(l1,state,L)
            else
                one_over_two = V_on_2/V_on
                d2 = rand()
                if d2 > one_over_two
                    Lift_1(l1,state,L)
                else
                    Lift_2(l2,state,L)
                end
            end
        elseif d1 > on_off_over_diffuse
            if V_off_1 == 0
                Detach(protein_position_2,state,l2)
            elseif V_off_2 == 0
                Detach(protein_position_1,state,l1)
            else
                one_over_two = V_off_2/V_off
                d2 = rand()
                if d2 > one_over_two
                    Detach(protein_position_1,state,l1)
                else 
                    Detach(protein_position_2,state,l2)
                end
            end
        else
            if V_d_1 == 0
                Diffuse_2(protein_position_2,state,L,l2)
            elseif V_d_2 == 0
                Diffuse_1(protein_position_1,state,L,l1)
            else
                one_over_two = V_d_2/V_d
                d2 = rand()
                if d2 > one_over_two
                    Diffuse_1(protein_position_1,state,L,l1)
                else
                    Diffuse_2(protein_position_2,state,L,l2)
                end
            end
        end
        Δt = Δt + δt
        t = t + δt
        push!(Time,t)
        push!(Length,sum([abs(i) for i in state])/(L))
        push!(state_1,sum([1 for i in state if i==1])/L)
        push!(state_2,sum([1 for i in state if i==-1])/L)
        for j = 1:(l1+l2)
            push!(Gaps[j],gaps_count(state,j))
        end
    end
end
export Simulate
end
