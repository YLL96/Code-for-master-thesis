M = 1000
base = 100
spec = 2000 / base
#function solve_dc_bi(thermalLimitscale, generators, buses, lines)
m = Model(with_optimizer(CPLEX.Optimizer))
@variable(m, plf[1:numlines])
@variable(m, vl[1:numlines], Bin)
@variable(m, δ[1:numbuses])
@variable(m, psd[1:numbuses])
@variable(m, pg[1:numgenerators])
@variable(m, g[1:numgenerators], Bin)
@variable(m, μ[1:numlines])
@variable(m, λ[1:numbuses])
@variable(m, θ_low[1:numgenerators] >= 0)
@variable(m, θ_up[1:numgenerators] >= 0)
@variable(m, ω_low[1:numlines] >= 0)
@variable(m, ω_up[1:numlines] >= 0)
@variable(m, α_low[1:numbuses] >= 0)
@variable(m, α_up[1:numbuses] >= 0)
@variable(m, θlow[1:numgenerators], Bin)
@variable(m, θup[1:numgenerators], Bin)
@variable(m, ωlow[1:numlines], Bin)
@variable(m, ωup[1:numlines], Bin)
@variable(m, αlow[1:numbuses], Bin)
@variable(m, αup[1:numbuses], Bin)

@objective(m, Max, sum(psd[i] for i=1:numbuses))   #objective function that maximizes generation cost
#@objective(m, Min, sum((1 - vl[i]) for i = 1:numlines))
#upper level problem
@constraint(m, sum((1-vl[i]) for i=1:numlines)<=1)
#@constraint(m, 10<=sum(psd[i] for i = 1:numbuses))

#lower level KKT


for i = 1:numbuses
    @constraint(m,sum(pg[j] for j in 1:numgenerators if genid[j] == i) + psd[i] ==sum(A[i, l] * plf[l] for l = 1:numlines) + Pd[i])
    @constraint(m, 0 <= psd[i] <= Pd[i])
    @constraint(m, sum(βl[i] * A[i, l] * μ[l] for l = 1:numlines) == 0)
    @constraint(m, 1 - λ[i] - α_low[i] + α_up[i] == 0)
    @constraint(m, α_low[i] <= M * αlow[i])
    @constraint(m, psd[i] <= M * (1 - αlow[i]))
    @constraint(m, α_up[i] <= M * αup[i])
    @constraint(m, Pd[i] - psd[i] <= M * (1 - αup[i]))

end


for j = 1:numlines
    @constraint(m, plf[j] == βl[j] * sum(A[n, j] * δ[n] for n = 1:numbuses))
    @constraint(m, 0 <= plf[j] + Plmax[j] * vl[j])#thermal flow limit
    @constraint(m, 0 <= Plmax[j] * vl[j] - plf[j])
    @constraint(m,sum(A[n, j] * λ[n] for n = 1:numbuses) - μ[j] - ω_low[j] + ω_up[j] == 0)
    @constraint(m, ω_low[j] <= M * ωlow[j])
    @constraint(m, plf[j] + Plmax[j] * vl[j] <= M * (1 - ωlow[j]))
    @constraint(m, ω_up[j] <= M * ωup[j])
    @constraint(m, Plmax[j] * vl[j] - plf[j] <= M * (1 - ωup[j]))

end

for i = 1:numgenerators
    @constraint(m,  Pgmin[i]<= pg[i])
    @constraint(m, pg[i] <= Pgmax[i])
    #关于genids的大问题
    @constraint(m,-sum(λ[n] for n = 1:numbuses if n==genid[i]) - θ_low[i] + θ_up[i] ==0)
    @constraint(m, θ_low[i] <= M * θlow[i])
    @constraint(m, pg[i] - Pgmin[i] <= M * (1 - θlow[i]))
    @constraint(m, θ_up[i] <= M * θup[i])
    @constraint(m, Pgmax[i] - pg[i] <= M * (1 - θup[i]))

end


optimize!(m)

optimal_objective = objective_value(m)
