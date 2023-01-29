using Distributed
addprocs(6)
N = 6
@everywhere begin
    using ProgressMeter
end
cmds = []
for i in 11:11+N
    cmd = `julia fret_stats_refined_fast_implementation.jl $i`
    push!(cmds,cmd)
end
@showprogress 1 "sampling" pmap(run,cmds)