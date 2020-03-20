m =Model(with_optimizer(Ipopt.Optimizer))
@variable(m, sij[1:numlines])
@variable(m, θbar[1:numbuses])
@variable(m, -θᵘ <= θij[1:numlines] <= θᵘ)
@variable(m, pgbar[1:numgenerators])
@variable(m, αg[1:numgenerators])
@variable(m, δ[1:numbuses])
@objective(m, Min, sum(c0[i]+c1[i]*pgbar[i]+c2[i]*(pgbar[i]^2+αg[i]^2*total_var) for i=1:numgenerators))


for i=2:numbuses

	@constraint(m, sum(B[i,j]δ[j] for j=2:numbuses)-sum(αg[j] for j=1:numgenerators if genid[j]==i)==0.)
	@constraint(m, sum(B[i,j]θbar[j] for j=2:numbuses)==fn[i]+sum(pgbar[j] for j=1:numgenerators if genid[j]==i)-Pd[i])

end

@constraint(m, sum(B[refbus,j]θbar[j] for j=1:numbuses)==-Pd[refbus])
@constraint(m, sum(αg[k] for k=1:numgenerators)==1)
@constraint(m,δ[refbus]==0.0 )
@constraint(m,θbar[refbus]==0.0)
#本不应该有的多余条件
#@constraint(m, sum(pgbar[i] for i=1:numgenerators)+sum(fn[i] for i in 1:numbuses)==sum(Pd[i] for i=1:numbuses))
for i=1:numlines

	@constraint(m, θij[i]==θbar[head[i]]-θbar[tail[i]])
	#@constraint(m, βl[i]*θij[i]<=Plmax[i])
	#@constraint(m, -Plmax[i]<=βl[i]*θij[i])
	@constraint(m, βl[i]*θij[i]+βl[i]*2.326*sij[i]^2<=Plmax[i])
	@constraint(m, -Plmax[i]<=βl[i]*θij[i]-βl[i]*2.326*sij[i]^2)
	@constraint(m, (sum(farm[k].σ^2*(𝚷[head[i],farm[k].bus]-𝚷[tail[i],farm[k].bus]-δ[head[i]]+δ[tail[i]])^2 for k=1:numfarms))<=sij[i]^2)
	@constraint(m, 0<=sij[i])

end

for i=1:numgenerators
	@constraint(m, -αg[i]<=0)
    @constraint(m, -pgbar[i]<=0)
	@constraint(m, pgbar[i]<=Pgmax[i])
	@constraint(m, Pgmin[i]<=pgbar[i])
	#@constraint(m, -pgbar[i]<=-Pgmin[i]-3.719*total_var*αg[i])
	#@constraint(m, pgbar[i]<=Pgmax[i]-3.719*total_var*αg[i])
end

optimize!(m)

optimal_objective = objective_value(m)
pgbar_opt = value.(pgbar)
αg_opt = value.(αg)
θbar_opt=value.(θbar)
δ_opt=value.(δ)
println("*************************************************************************")
println("The optimal cost is ", optimal_objective)
println("*************************************************************************")
