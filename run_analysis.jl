using Distributed
addprocs(10)
# exp_label="wt_15mM_salt"
# simu_label="with_diffusion_length"
# Ls = [1000]
# for L in Ls
#     run(`julia evaluate.jl $(exp_label) $(simu_label)_$(L)`)
# end
# for L in Ls
#     run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)_$(L) `)
# end
# exp_label="wt_15mM_salt"
# simu_label="length"
# Ls = [300,800,1000,5000]
# for L in Ls
#     run(`julia evaluate.jl $(exp_label) $(simu_label)_$(L)`)
# end
# for L in Ls
#     run(`julia exact_gap_analysis.jl $(exp_label) $(simu_label)_$(L) `)
# end

exp_labels = ["wt_15mM_salt","wt_150mM_salt"]
simu_label = "init"
Ls = [1000]
cmds = Array{Cmd,1}()
for exp_label in exp_labels
    for L in Ls
        push!(cmds,`julia evaluate.jl $(exp_label) $(simu_label) `)
        # push!(cmds,`julia exact_gap_analysis.jl $(exp_label) $(simu_label) `)
    end
end
pmap(run,cmds)