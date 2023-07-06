using DifferentialEquations, Distributions
# Define the base ODEProblem
f(u,p,t) = 0.98u
u0 = 1.0 # some initial value
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
# Define the prob_func that modifies the initial value
function prob_func(prob,i,repeat)
  remake(prob,u0=rand(Normal(1.0,0.1))) # sample a new u0 from N(1, 0.1)
end
# Define the EnsembleProblem with 100 trajectories
ens_prob = EnsembleProblem(prob, prob_func=prob_func)
# Solve the EnsembleProblem
sim = solve(ens_prob,Tsit5(),EnsembleThreads(),trajectories=100)
