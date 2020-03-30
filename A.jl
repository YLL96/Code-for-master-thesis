using JuMP
using CPLEX
using Distributions
using PowerModels
root="/Users/yuyangli/Desktop/matpower6.0/"
#casename="case24_ieee_rts.m"
casename="case118.m"
path=root*casename
network_data=PowerModels.parse_file(path)
result=run_dc_opf(path, with_optimizer(CPLEX.Optimizer))
numbuses=length(collect(keys(network_data["bus"])))
numgenerators=length(collect(keys(network_data["gen"])))
numlines=length(collect(keys(network_data["branch"])))
numloads=length(collect(keys(network_data["load"])))
println(numgenerators, " generators ", numbuses,  " buses ", numlines,  " lines ")

loadscale=1.
n=numbuses
l=numlines
A= zeros(n,l)
for l in 1:numlines
    count=string(l)
    i = network_data["branch"][count]["f_bus"]
    j = network_data["branch"][count]["t_bus"]
    A[i,l] = 1
    A[j,l] = -1
end

#Maximum and Minimum generator output& genID
Pgmax=zeros(numgenerators)
Pgmin=zeros(numgenerators)
genid=[]
c0=zeros(numgenerators)
c1=zeros(numgenerators)
c2=zeros(numgenerators)
for g in 1:numgenerators
    count=string(g)
    Pgmax[g]=network_data["gen"][count]["pmax"]
    Pgmin[g]=network_data["gen"][count]["pmin"]
    if isempty(network_data["gen"][count]["cost"])
        c2[g]=0
        c1[g]=0
        c0[g]=0
    elseif length(network_data["gen"][count]["cost"])==2
        c2[g]=0
        c1[g]=network_data["gen"][count]["cost"][1]
        c0[g]=network_data["gen"][count]["cost"][2]
    else
        c2[g]=network_data["gen"][count]["cost"][1]
        c1[g]=network_data["gen"][count]["cost"][2]
        c0[g]=network_data["gen"][count]["cost"][3]
    end
    push!(genid,network_data["gen"][count]["gen_bus"])
end

#Thermal Rating of lines&susceptance
Plmax=zeros(numlines)
βl=zeros(numlines)
head=[]
tail=[]
#case118
lmax= [ 175 175	500	175	175	175	500	500	500	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	500	500	500	175	175	500	175	500	175	175	140	175	175	175	175	175	175	175	175	500	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	500	500	500	500	500	500	500	175	175	500	175	500	175	175	500	500	175	175	175	175	175	175	175	500	175	175	175	175	175	175	500	500	175	500	500	200	200	175	175	175	500	500	175	175	500	500	500	175	500	500	175	175	175	175	175	175	175	175	175	175	200	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	500	175	175	175 ]
for l in 1:numlines
    count=string(l)
    if haskey(network_data["branch"][count],["rate_a"])==true
        Plmax[l]=network_data["branch"][count]["rate_a"]
    else
        Plmax[l]=lmax[l]
    end
    βl[l]=1/network_data["branch"][count]["br_x"]
    push!(head, network_data["branch"][count]["f_bus"])
    push!(tail, network_data["branch"][count]["t_bus"])
end
refbus=1

#node that send or take power flow from/to node i
take=Vector{Int}[]
give=Vector{Int}[]

for i in 1:numbuses
    temp_take=[]
    temp_give=[]
    for l in 1:numlines
        if i==head[l]
            push!(temp_take,tail[l])
        end
        if i==tail[l]
            push!(temp_give,head[l])
        end
    end
    push!(take, temp_take)
    push!(give, temp_give)
end

#Demand at each bus
Pd=zeros(numbuses)
for n in 1:numloads
    count=string(n)
    place=network_data["load"][count]["load_bus"]
    Pd[place]=network_data["load"][count]["pd"]
end

#B matrix
n = numbuses
B = zeros(n,n)
for l in 1:numlines
    count=string(l)
    i = network_data["branch"][count]["t_bus"]
    j = network_data["branch"][count]["f_bus"]
    B[i,j] = -βl[l]
    B[j,i] = -βl[l]
    B[i,i] += βl[l]
    B[j,j] += βl[l]
end

mutable struct Farm
    μ::Float64
    σ::Float64
    bus::Int
end
#Data for windfarm
#case39
wp = 0.0
factor_σ = wp*0.
farm=Farm[]
push!(farm, Farm(113.0/100*wp, 	factor_σ *11.3/100, 24))
push!(farm, Farm(84.0/100*wp,   factor_σ *8.4/100, 	26))
push!(farm, Farm(250.0/100*wp, 	factor_σ *25.0/100, 38))
push!(farm, Farm(118.0/100*wp,  factor_σ *11.8/100, 43))
#push!(farm, Farm(76.0/100*wp,  	factor_σ *7.6/100, 	49))
#push!(farm, Farm(72.0/100*wp,  	factor_σ *7.2/100, 	53))

total_var=sum(f.σ^2 for f in farm)
numfarms=length(farm)
fn=zeros(numbuses)
stdn=zeros(numbuses)
for n in 1:numbuses
    for w in 1:numfarms
        if n==farm[w].bus
            fn[n]=farm[w].μ
            stdn[n]=farm[w].σ
        end
    end
end
@assert size(B,1) == size(B,2)
n = size(B,1)
Bhat=B[1:n .!= refbus, 1:n .!= refbus]
Binv=inv(Bhat)
a=zeros(numbuses-1)
b=zeros((1,numbuses))
c=hcat(a,Binv)
𝚷=vcat(b,c)

theta_u= 15
θᵘ = deg2rad(float(theta_u))
rate=sum(fn)/sum(Pd)
#Check if Bθ=0
θr=[]
for i=1:numbuses
    count=string(i)
    push!(θr, result["solution"]["bus"][count]["va"])
end
sum(B*θr)
