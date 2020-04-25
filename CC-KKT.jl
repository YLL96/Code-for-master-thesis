include("A.jl")
#M=100000000
M=1000
#M=1000
percent=0.4
shed=sum(Pd)*percent
m =Model(CPLEX.Optimizer)
#@variable(m, sij[1:numlines])
@variable(m, θbar[1:numbuses])
@variable(m, pgbar[1:numgenerators])
@variable(m, αg[1:numgenerators])
@variable(m, δ[1:numbuses])
@variable(m, psd[1:numbuses])

@variable(m, vl[1:numlines], Bin)

@variable(m, μn[1:numbuses-1])
@variable(m, λn[1:numbuses-1])
@variable(m, y1)
@variable(m, y2)

@variable(m, y5)
@variable(m, y6)

@variable(m, e_l[1:numbuses]>=0)
@variable(m, e_u[1:numbuses]>=0)
@variable(m, ω1[1:numgenerators]>=0)
@variable(m, ω2[1:numgenerators]>=0)
@variable(m, m_l[1:numgenerators]>=0)
@variable(m, m_u[1:numgenerators]>=0)
#@variable(m, κ_l[1:numlines]>=0)
#@variable(m, κ_u[1:numlines]>=0)

@variable(m, c_l[1:numlines]>=0)
@variable(m, c_u[1:numlines]>=0)
#@variable(m, g_l[1:numlines]>=0)
#@variable(m, g_u[1:numlines]>=0)

@variable(m, w_e_l[1:numbuses], Bin)
@variable(m, w_e_u[1:numbuses], Bin)
@variable(m, w_m_l[1:numgenerators], Bin)
@variable(m, w_m_u[1:numgenerators], Bin)
#@variable(m, w_κ_l[1:numlines], Bin)
#@variable(m, w_κ_u[1:numlines], Bin)

@variable(m, w_c_l[1:numlines], Bin)
@variable(m, w_c_u[1:numlines], Bin)
#@variable(m, w_g_l[1:numlines], Bin)
#@variable(m, w_g_u[1:numlines], Bin)

@variable(m, w1[1:numgenerators], Bin)
@variable(m, w2[1:numgenerators], Bin)

#Minimum vulnerability analysis
@objective(m, Min, sum((1 - vl[i]) for i = 1:numlines))
@constraint(m, shed<=sum(psd[i] for i = 1:numbuses))
#Maximum vulnerability analysis
#@objective(m, Max, sum(psd[i] for i=1:numbuses))
#@constraint(m, sum((1 - vl[i]) for i = 1:numlines)<=1)
#lower level KKT
for i= 1:numbuses-1

  @constraint(m, sum(B[i,j]δ[j] for j=1:numbuses-1)-sum(αg[j] for j=1:numgenerators if genid[j]==i)==0.)
  @constraint(m, sum(B[i,j]θbar[j] for j=1:numbuses-1)==fn[i]+psd[i]+sum(pgbar[j] for j=1:numgenerators if genid[j]==i)-Pd[i])

  @constraint(m, psd[i]>=0)
  @constraint(m, e_l[i]<=M*w_e_l[i])
  @constraint(m, psd[i]<=M*(1-w_e_l[i]))

  @constraint(m, psd[i]<=Pd[i])
  @constraint(m, e_u[i]<=M*w_e_u[i])
  @constraint(m, -psd[i]+Pd[i]<=M*(1-w_e_u[i]))
  #Δpd
  @constraint(m, 1+μn[i]-y1-e_l[i]+e_u[i]==0)
  #δ
  @constraint(m, sum(-B[m,i]λn[m] for m=1:numbuses-1)==0)
  #θ
  #@constraint(m, sum(-B[m,i]μn[m] for m=1:numbuses-1)+sum(-κ_l[l]βl[l]+κ_u[l]βl[l]-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if head[l]==i)-sum(-κ_l[l]βl[l]+κ_u[l]βl[l]-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if tail[l]==i)==0)
  @constraint(m, sum(-B[m,i]μn[m] for m=1:numbuses-1)+sum(-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if head[l]==i)-sum(-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if tail[l]==i)==0)


end

#System wide balance
@constraint(m, sum(pgbar[i] for i=1:numgenerators)+sum(fn[i] for i in 1:numbuses)==sum(Pd[i]-psd[i] for i=1:numbuses))
@constraint(m, sum(αg[k] for k=1:numgenerators)==1.)
@constraint(m,δ[refbus]==0.0 )
@constraint(m,θbar[refbus]==0.0)

@constraint(m, 0<=psd[refbus]<=Pd[refbus])
@constraint(m, e_l[refbus]<=M*w_e_l[refbus])
@constraint(m, psd[refbus]<=M*(1-w_e_l[refbus]))
@constraint(m, e_u[refbus]<=M*w_e_u[refbus])
@constraint(m, -psd[refbus]+Pd[refbus]<=M*(1-w_e_u[refbus]))
#Δpd
@constraint(m, 1-y1-e_l[refbus]+e_u[refbus]==0)

#θ
#@constraint(m, -y6+sum(-κ_l[l]βl[l]+κ_u[l]βl[l]-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if head[l]==refbus)-sum(-κ_l[l]βl[l]+κ_u[l]βl[l]-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if tail[l]==refbus)==0)
@constraint(m, -y6+sum(-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if head[l]==refbus)-sum(-c_l[l]βl[l]+c_u[l]βl[l] for l=1:numlines if tail[l]==refbus)==0)

#δ
@constraint(m, -y5==0)





for i=1:numlines
#=
  @constraint(m, βl[i]*(θbar[head[i]]-θbar[tail[i]])<=Plmax[i]*vl[i])
  @constraint(m, κ_l[i]<=M*w_κ_l[i])
  @constraint(m, vl[i]*Plmax[i]-βl[i]*(θbar[head[i]]-θbar[tail[i]])<=M*(1-w_κ_l[i]))

  @constraint(m, -Plmax[i]*vl[i]<=βl[i]*(θbar[head[i]]-θbar[tail[i]]))
  @constraint(m, κ_u[i]<=M*w_κ_u[i])
  @constraint(m, βl[i]*(θbar[head[i]]-θbar[tail[i]])+vl[i]*Plmax[i]<=M*(1-w_κ_u[i]))
=#
  @constraint(m, -Plmax[i]*vl[i]<=βl[i]*(θbar[head[i]]-θbar[tail[i]])-βl[i]*su*vl[i]*varlimit[i])
  @constraint(m, c_l[i]<=M*w_c_l[i])
  @constraint(m, -su*βl[i]*vl[i]*varlimit[i]+Plmax[i]*vl[i]+βl[i]*(θbar[head[i]]-θbar[tail[i]])<=M*(1-w_c_l[i]))

  @constraint(m, βl[i]*(θbar[head[i]]-θbar[tail[i]])+βl[i]*su*varlimit[i]*vl[i]<=Plmax[i]*vl[i])
  @constraint(m, c_u[i]<=M*w_c_u[i])
  @constraint(m, Plmax[i]*vl[i]-βl[i]*(θbar[head[i]]-θbar[tail[i]])-su*βl[i]*vl[i]*varlimit[i]<=M*(1-w_c_u[i]))
#=
  @constraint(m, varlimit[i]*vl[i]<=sij[i])
  @constraint(m, g_l[i]<=M*w_g_l[i])
  @constraint(m, -vl[i]*varlimit[i]+sij[i]<=M*(1-w_g_l[i]))

  @constraint(m, sij[i]<=M*vl[i])
  @constraint(m, g_u[i]<=M*w_g_u[i])
  @constraint(m, M*vl[i]-sij[i]<=M*(1-w_g_u[i]))

  #sij/hij
  @constraint(m, su*c_l[i]βl[i]+su*c_u[i]βl[i]-g_l[i]+g_u[i]==0)
  #@constraint(m, su*c_l[i]βl[i]+su*c_u[i]βl[i]-g_l[i]==0)
=#

end

for i=1:numgenerators

  @constraint(m, αg[i]>=0)
  @constraint(m, ω1[i]<=M*w1[i])
  @constraint(m, αg[i]<=M*(1-w1[i]))

  @constraint(m, pgbar[i]>=0)
  @constraint(m, ω2[i]<=M*w2[i])
  @constraint(m, pgbar[i]<=M*(1-w2[i]))

  @constraint(m, -pgbar[i]<=-Pgmin[i]-su*total_var*αg[i])
  @constraint(m, m_l[i]<=M*w_m_l[i])
  @constraint(m, pgbar[i]-Pgmin[i]-su*total_var*αg[i]<=M*(1-w_m_l[i]))

  @constraint(m, pgbar[i]<=Pgmax[i]-su*total_var*αg[i])
  @constraint(m, m_u[i]<=M*w_m_u[i])
  @constraint(m, -pgbar[i]+Pgmax[i]-su*total_var*αg[i]<=M*(1-w_m_u[i]))
  #pgbar
  @constraint(m, -y1+sum(μn[n] for n=1:numbuses-1 if n==genid[i])-ω2[i]-m_l[i]+m_u[i]==0)
  #αg
  @constraint(m, sum(λn[n] for n=1:numbuses-1 if genid[i]==n)-y2-ω1[i]+su*total_var*m_l[i]+su*total_var*m_u[i]==0)


end

optimize!(m)

optimal_objective = objective_value(m)
pgbar_opt = value.(pgbar)
αg_opt = value.(αg)
θbar_opt=value.(θbar)
δ_opt=value.(δ)
vl_opt=value.(vl)
psd_opt=value.(psd)
y5_opt=value.(y5)
#sij_opt=value.(sij)
println("*************************************************************************")
println("The optimal cost is ", optimal_objective)
println("*************************************************************************")
pl=[]
for l=1:numlines
  temp_pl=βl[l]*(θbar_opt[head[l]]-θbar_opt[tail[l]])
  push!(pl, temp_pl)
end
