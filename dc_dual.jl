m =Model(with_optimizer(CPLEX.Optimizer))
#@variable(m, plf[1:numlines])
@variable(m, vl[1:numlines], Bin)
@variable(m, μ[1:numlines])
@variable(m, λ[1:numbuses])
@variable(m, obj)
@variable(m, θ_low[1:numgenerators]>=0)
@variable(m, θ_up[1:numgenerators]>=0)
@variable(m, ω_low[1:numlines]>=0)
@variable(m, ω_up[1:numlines]>=0)
@variable(m, α_low[1:numbuses]>=0)
@variable(m, α_up[1:numbuses]>=0)

@variable(m, tl[1:numlines])
@variable(m, hl[1:numlines])

@variable(m, t_low[1:numlines])
@variable(m, t_up[1:numlines])
@variable(m, h_low[1:numlines])
@variable(m, h_up[1:numlines])


@objective(m, Min, sum((1 - vl[i]) for i = 1:numlines))#objective function that maximizes generation cost
@constraint(m, obj==sum((λ[n]-α_up[n])Pd[n] for n=1:numbuses)-sum((ω_low[l]+ω_up[l])Plmax[l] for l=1:numlines)+sum(θ_low[g]Pgmin[g]-θ_up[g]Pgmax[g] for g=1:numgenerators))
@constraint(m, 6 <=obj)

for i=1:numbuses
    @constraint(m, sum(βl[i]*A[i,l]*tl[l] for l = 1:numlines)==0)
    @constraint(m, 1-λ[i]-α_low[i]+α_up[i]==0)
end
for i=1:numgenerators
    @constraint(m, -sum(λ[n] for n = 1:numbuses if n==genid[i])-θ_low[i]+θ_up[i]==0)
end
for i=1:numlines
    @constraint(m, sum(A[n,i]λ[n] for n=1:numbuses)-μ[i]-ω_low[i]+ω_up[i]==0)
    #这里有问题

    @constraint(m, tl[i]==μ[i]-hl[i])
    @constraint(m, -1*vl[i]<=tl[i])
    @constraint(m, tl[i]<=1*vl[i])
    @constraint(m, -1*(1-vl[i])<=hl[i])
    @constraint(m, hl[i]<=1*(1-vl[i]))

end
optimize!(m)
optimal_objective = objective_value(m)
vl=JuMP.value.(vl)
