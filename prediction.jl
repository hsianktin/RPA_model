# prediction
using DataFrames
using Distributed
using CSV
addprocs(12)

@everywhere begin
    using ProgressMeter
end

cmds = Array{Cmd,1}()
L = 1000

exp_label = "wt_15mM_salt"
simu_label = "prediction"
paras = [1e-5,1e-3, 1e-4,1e-3,2,3.2]
folds = [50,100]
N = 100
k_on,k_off,k_open,k_close = paras[1:4]
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"

for fold in folds
    cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label) $gaps_type`
    run(cmd)
end


L = 1000

exp_label = "wt_150mM_salt"
simu_label = "prediction"
paras = [1e-5,1e-2, 1e-4,1e-2,2,3.2]
folds = [50,100]
N = 100
k_on,k_off,k_open,k_close = paras[1:4]
T1 = 1800.0
T2 = 600.0
gaps_type = "exact"
for fold in folds
    cmd=`julia simu_base.jl $k_on $k_off $k_open $k_close $fold $L $T1 $T2 $N $(exp_label) $(simu_label) $gaps_type`
    run(cmd)
end