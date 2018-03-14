# Let us first make the Convex and MultiConvex modules available
using Convex
using ECOS
using Gurobi
using SCS
using Ipopt
using PyPlot
using HDF5
using JuMP

REUSED_INPUT = true
FLOT_ALL_FIGS = true

F= 5*10^4
S= 220
B= 1000
L= 100*1024*8  #File size
No_MNO=3
max_iters = 100
ALG_MODE = 4 #1, 2, 3, 4
w_tradeoff= 1.

input_filename = "weights$No_MNO.h5"
if (REUSED_INPUT)
   h5open(input_filename, "r") do file
      global lamb_n = read(file, "lamb_n")
      global w_ulti = read(file, "w_ulti")
      global w_lamb_hit_n = read(file, "w_lamb_hit_n")
   end
else
  include("Setting$No_MNO.jl")
  h5open(input_filename, "w") do file
      write(file, "lamb_n", lamb_n)
      write(file, "w_ulti", w_ulti)
      write(file, "w_lamb_hit_n", w_lamb_hit_n)
  end
end

include("Plot_Fig_general.jl")

println("W_util:",w_ulti)
println("lamb_n",lamb_n)
println("W_lamb_hit_n:",w_lamb_hit_n)

n=No_MNO
epsilon = 1e-3

alphas = zeros(n, max_iters)
betas  = zeros(n, max_iters)
gammas = zeros(n, max_iters)
total_cost = zeros(max_iters)
langragians = zeros(max_iters)
s_fes = zeros(n,max_iters)
dual_vars1 = zeros(max_iters)
dual_vars2 = zeros(max_iters)
dual_vars3 = zeros(max_iters)
max_w    = 1e6
w_factor = 1.2
w_init   = 0.1
# w_init   = 0.5

# rho = 0.5
# rho1 = 15.
rho = 0.5
rho1 = 24.
Jacobian_step = 1.
if ALG_MODE > 2
  delta = 2.
  rho = (rho1*(No_MNO/(2-Jacobian_step)-1) + delta)
end

function beta_obj(beta,lamb_miss, k)
  obj  =  vecdot(w_ulti,beta./alphas[:,k]) + rho/2*sumsquares(beta - betas[:,k]) +
  # w_tradeoff*(-No_MNO + sum((gammas[:,k]*B).*invpos(gammas[:,k]*B-lamb_miss)) )
  w_tradeoff*(-No_MNO + sum( max(0,(gammas[:,k]*B).*invpos(gammas[:,k]*B-lamb_miss)))  )
  # w_tradeoff*sum( max(lamb_miss,3*lamb_miss -2/3*gammas[:,k]*B,10*lamb_miss - 16/3*gammas[:,k]*B,
  # 70*lamb_miss - 178/3*gammas[:,k]*B, 500*lamb_miss - 1468/3*gammas[:,k]*B,5000*lamb_miss - 16318/3*gammas[:,k]*B ) )

  return obj
end

function beta_i_obj(beta,lamb_miss, i, k)
  sub_obj = beta
  for idx = 1:n
      if idx< i
        # sub_obj += betas[idx,k+1] # Gauss-Seidel
        sub_obj += betas[idx,k]     # Jacobian
      elseif idx > i
        sub_obj += betas[idx,k]
      end
  end

  obj  =  w_ulti[i]*beta/alphas[i,k] + rho/2*square(beta - betas[i,k]) + rho1/2*square(sub_obj - 1 - dual_vars1[k]/rho1) +
  w_tradeoff*(-1 + max(0,gammas[i,k]*B*invpos(gammas[i,k]*B-lamb_miss)) )
  # max( lamb_miss,3*lamb_miss -2/3*gammas[i,k]*B,10*lamb_miss - 16/3*gammas[i,k]*B,70*lamb_miss - 178/3*gammas[i,k]*B,
  # 500*lamb_miss - 1468/3*gammas[i,k]*B,5000*lamb_miss - 16318/3*gammas[i,k]*B )

  return obj
end


function alpha_gamma_obj(alpha,gamma,lamb_miss, k)
  obj  = vecdot(w_ulti,betas[:,k+1].*invpos(alpha)) + rho/2*sumsquares(alpha - alphas[:,k]) + rho/2*sumsquares(gamma - gammas[:,k]) +
  # w_tradeoff*sum(lamb_miss.*invpos(gamma*B-lamb_miss))
  w_tradeoff*sum( max(0,lamb_miss.*invpos(gamma*B-lamb_miss)) )
  # w_tradeoff*sum( max( lamb_miss,3*lamb_miss -2/3*gamma*B,10*lamb_miss - 16/3*gamma*B,70*lamb_miss - 178/3*gamma*B,
  # 500*lamb_miss - 1468/3*gamma*B,5000*lamb_miss - 16318/3*gamma*B ) )

  return obj
end

function alpha_gamma_i_obj(alpha,gamma,lamb_miss,i, k)
  sub_obj2 = alpha
  sub_obj3 = gamma
  for idx = 1:n
      if idx< i
        # sub_obj2 += alphas[idx,k+1]  # Gauss-Seidel
        # sub_obj3 += gammas[idx,k+1]  # Gauss-Seidel
        sub_obj2 += alphas[idx,k]  # Jacobian
        sub_obj3 += gammas[idx,k]  # Jacobian
      elseif idx > i
        sub_obj2 += alphas[idx,k]
        sub_obj3 += gammas[idx,k]
      end
  end

  obj  = w_ulti[i]*betas[i,k+1]*invpos(alpha)+ rho/2*((alpha - alphas[i,k])^2)+ rho/2*((gamma - gammas[i,k])^2) +
  rho1/2*(sub_obj2 - 1 - dual_vars2[k]/rho1)^2 + rho1/2*(sub_obj3 - 1 - dual_vars3[k]/rho1)^2 +
  # w_tradeoff*(lamb_miss*invpos(gamma*B-lamb_miss))
  w_tradeoff*max(0,lamb_miss*invpos(gamma*B-lamb_miss) )
  # max( lamb_miss,3*lamb_miss -2/3*gamma*B,10*lamb_miss - 16/3*gamma*B,70*lamb_miss - 178/3*gamma*B,
  # 500*lamb_miss - 1468/3*gamma*B, 5000*lamb_miss - 16318/3*gamma*B )

  return obj
end

function beta_i_update0(i,k)
  beta = Convex.Variable()
  lamb_miss = lamb_n[i] - w_lamb_hit_n[i] * beta
  problem = minimize(beta_i_obj(beta,lamb_miss, i, k))
  problem.constraints += [lamb_miss>=0, beta>=0, w_ulti[i]*beta  - (1 - epsilon)*alphas[i,k] <=0,
                          lamb_miss < gammas[i,k]*B]

  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_i_update0(i,k)
  alpha = Convex.Variable()
  gamma = Convex.Variable()
  lamb_miss = lamb_n[i] - w_lamb_hit_n[i] *betas[i,k+1]
  problem = minimize(alpha_gamma_i_obj(alpha,gamma,lamb_miss,i, k))
  problem.constraints += [alpha > 0, gamma >= 0, w_ulti[i]*betas[i,k+1] - (1 - epsilon)*alpha <=0,
                          lamb_miss < gamma*B]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("alphagamma-subprob: ",problem.status)
  # println("obj:", problem.optval)
  # if i==1
  #   println("obj1:", w_ulti*betas[i,k+1]/(alpha.value) - gamma.value + rho*(alpha.value - alphas[i,k])^2+ rho*(gamma.value - gammas[i,k])^2
  #   + rho1*(alpha.value+ alphas[2,k] - 1 - dual_vars2[k])^2 + rho1*(gamma.value+ gammas[2,k]- 1 - dual_vars3[k])^2 )
  #   println("alpha deviation1:", rho1*(alpha.value+ alphas[2,k] - 1 - dual_vars2[k])^2 )
  #   println("gamma deviation1:", rho1*(gamma.value+ gammas[2,k]- 1 - dual_vars3[k])^2 )
  # elseif i==2
  #   println("obj2:", w_ulti*betas[i,k+1]/(alpha.value) - gamma.value + rho*(alpha.value - alphas[i,k])^2+ rho*(gamma.value - gammas[i,k])^2
  #   + rho1*(alpha.value+ alphas[1,k+1] - 1 - dual_vars2[k])^2 + rho1*(gamma.value+ gammas[1,k+1]- 1 - dual_vars3[k])^2 )
  #   println("alpha deviation2:", rho1*(alpha.value+ alphas[1,k+1] - 1 - dual_vars2[k])^2)
  #   println("gamma deviation2:", rho1*(gamma.value+ gammas[1,k+1]- 1 - dual_vars3[k])^2)
  # end
  return alpha.value, gamma.value, 0
end

function beta_i_update(i,k)
  beta = Convex.Variable()
  s = Convex.Variable()
  lamb_miss = lamb_n[i] - w_lamb_hit_n[i] * beta
  problem = minimize(beta_i_obj(beta,lamb_miss, i, k) + w*s)
  problem.constraints +=  [lamb_miss>=0,beta>=0, w_ulti[i]*beta  - (1 - epsilon)*alphas[i,k] <=s, s>=0,
                          lamb_miss < gammas[i,k]*B]

  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_i_update(i,k)
  alpha = Convex.Variable()
  gamma = Convex.Variable()
  s = Convex.Variable()
  lamb_miss = lamb_n[i] - w_lamb_hit_n[i] *betas[i,k+1]
  problem = minimize(alpha_gamma_i_obj(alpha,gamma,lamb_miss,i, k) + w*s )
  problem.constraints +=  [alpha > 0, gamma >= 0, w_ulti[i]*betas[i,k+1] - (1 - epsilon)*alpha <=s, s>=0,
                          lamb_miss < gamma*B]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  return alpha.value, gamma.value, s.value
end

function beta_update0(k)
  beta = Convex.Variable(n)
  lamb_miss = lamb_n - w_lamb_hit_n .*beta
  problem = minimize(beta_obj(beta,lamb_miss, k))
  problem.constraints +=  [lamb_miss>=0, sum(beta) == 1, beta>=0, w_ulti.*beta  - (1 - epsilon)*alphas[:,k] <= 0,
                          lamb_miss < gammas[:,k]*B]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)
  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  return beta.value
end

function alpha_gamma_update0(k)
  alpha = Convex.Variable(n)
  gamma = Convex.Variable(n)
  lamb_miss = lamb_n - w_lamb_hit_n .*betas[:,k+1]
  problem = minimize(alpha_gamma_obj(alpha,gamma,lamb_miss, k))
  problem.constraints +=  [sum(alpha) == 1,sum(gamma) == 1,alpha > 0,gamma >= 0, w_ulti.*betas[:,k+1] - (1 - epsilon)*alpha <=0,
                          lamb_miss < gamma*B]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("alphagamma-subprob: ",problem.status)
  # println("alpha:", alpha.value)
  # println("gamma:", gamma.value)
  return alpha.value, gamma.value, zeros(n)
end

function beta_update(k)
  beta = Convex.Variable(n)
  s = Convex.Variable(n)
  lamb_miss = lamb_n - w_lamb_hit_n .*beta
  problem = minimize(beta_obj(beta,lamb_miss, k)+ w*sum(s))
  problem.constraints +=  [lamb_miss>=0,sum(beta) == 1, beta>=0, w_ulti.*beta  - (1 - epsilon)*alphas[:,k] <= s, s>= 0,
                          lamb_miss < gammas[:,k]*B]
  # solve!(prob, GurobiSolver(),verbose=false)
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)
  # solve!(prob, SCSSolver(verbose=false),verbose=false)
  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_update(k)
  alpha = Convex.Variable(n)
  gamma = Convex.Variable(n)
  s = Convex.Variable(n)
  lamb_miss = lamb_n - w_lamb_hit_n .*betas[:,k+1]
  problem = minimize(alpha_gamma_obj(alpha,gamma,lamb_miss, k)+ w*sum(s))
  problem.constraints +=  [sum(alpha) == 1,sum(gamma) == 1,alpha > 0,gamma >= 0, w_ulti.*betas[:,k+1] - (1 - epsilon)*alpha <= s, s>= 0,
                          lamb_miss < gamma*B]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("alphagamma-subprob: ",problem.status)
  # println("alpha:", alpha.value)
  # println("gamma:", gamma.value)
  # println("s:", s.value)
  return alpha.value, gamma.value, s.value
end

function cal_total_cost(k)
  lamb_miss = max(0,lamb_n - w_lamb_hit_n .*betas[:,k])
  # backhaul_cost = 0
  # for i = 1:n
  #   backhaul = gammas[i,k]*B
  #   fraction = lamb_miss[i]/backhaul
  #   if fraction < 1/3
  #     backhaul_cost += lamb_miss[i]
  #   elseif fraction < 2/3
  #     backhaul_cost += 3*lamb_miss[i] - 2/3*backhaul
  #   elseif fraction < 9/10
  #     backhaul_cost += 10*lamb_miss[i] - 16/3*backhaul
  #   elseif fraction < 1
  #     backhaul_cost += 70*lamb_miss[i] - 178/3*backhaul
  #   elseif fraction < 11/10
  #     backhaul_cost += 500*lamb_miss[i] - 1468/3*backhaul
  #   else
  #     backhaul_cost += 5000*lamb_miss[i] - 16318/3*backhaul
  #   end
  # end
  backhaul_cost = sum(lamb_miss./(gammas[:,k]*B - lamb_miss) )
  for i=1:No_MNO
    if gammas[i,k]*B < lamb_miss[i]
      println("M/M/1 Violation")
    end
  end

  total_cost[k] =  vecdot(w_ulti,betas[:,k]./alphas[:,k]) + w_tradeoff*backhaul_cost

  if(k==max_iters)
    println("Test lamb_miss:", lamb_miss)
    println("Test utilization:",vecdot(w_ulti,betas[:,k]./alphas[:,k]))
    println("Test miss cost:", lamb_miss./(gammas[:,k]*B - lamb_miss) )
    # println("Test miss cost1:", max(lamb_miss,3*lamb_miss -2/3*gammas[:,k]*B,10*lamb_miss - 16/3*gammas[:,k]*B,
    # 70*lamb_miss - 178/3*gammas[:,k]*B, 500*lamb_miss - 1468/3*gammas[:,k]*B,5000*lamb_miss - 16318/3*gammas[:,k]*B ) )
  end
end

# function cal_langragian(k)
#   langragians[k] = sum(betas[:,k]/alphas[:,k]) + sum(betas[:,k]) - sum(gammas[:,k]) + rho*(norm(alphas[:,k] - alphas[:,k-1])^2)
#    + rho*(norm(gammas[:,k] - gammas[:,k-1])^2) + rho*(norm(betas[:,k] - betas[:,k-1])^2)
# end

function save_results()
  h5open(string(folder,"results$ALG_MODE.h5"), "w") do file
     write(file, "alphas", alphas)
     write(file, "betas" , betas)
     write(file, "gammas", gammas)
     write(file, "total_cost", total_cost)
  end
end

function main()
  # alphas[:,1] = 1/No_MNO *ones(No_MNO)
  # betas[:,1]  = 1/No_MNO *ones(No_MNO)
  # gammas[:,1] = 1/No_MNO *ones(No_MNO)
  alphas[:,1] = 0.2*ones(No_MNO)
  betas[:,1]  = 0.2*ones(No_MNO)
  gammas[:,1] = 0.2*ones(No_MNO)
  cal_total_cost(1)
  # global rho = 0.5
  global w = w_init


  for k=1:max_iters-1
    # rho = rho
    w = min(max_w,w *w_factor)
    if ALG_MODE == 1
      betas[:,k+1]  = beta_update0(k)
      alphas[:,k+1], gammas[:,k+1], s_fes[:,k+1] = alpha_gamma_update0(k)
    elseif ALG_MODE == 2
      betas[:,k+1]  = beta_update(k)
      alphas[:,k+1], gammas[:,k+1], s_fes[:,k+1] = alpha_gamma_update(k)
    elseif ALG_MODE > 2
      for i =1:n
        if ALG_MODE == 3
          betas[i,k+1]  = beta_i_update0(i,k)
          alphas[i,k+1], gammas[i,k+1], s_fes[i,k+1] = alpha_gamma_i_update0(i,k)
        elseif ALG_MODE == 4
          betas[i,k+1]  = beta_i_update(i,k)
          alphas[i,k+1], gammas[i,k+1], s_fes[i,k+1] = alpha_gamma_i_update(i,k)
        end
      end
      dual_vars1[k+1] = dual_vars1[k] - rho1*Jacobian_step*(sum(betas[:,k+1])  - 1)
      dual_vars2[k+1] = dual_vars2[k] - rho1*Jacobian_step*(sum(alphas[:,k+1]) - 1)
      dual_vars3[k+1] = dual_vars3[k] - rho1*Jacobian_step*(sum(gammas[:,k+1]) - 1)
    end
    cal_total_cost(k+1)
  end

  plt_figures()
  save_results()

  println("alpha",alphas[:,max_iters-1:end])
  println("beta",betas[:,max_iters-1:end])
  println("gamma",gammas[:,max_iters-1:end])
  println("total_cost",total_cost[max_iters-1:end])
  println("s_fes:",s_fes[:,max_iters-1:end])
  println("constraint1:",w_ulti[1]*betas[1,max_iters-1:end]/alphas[1,max_iters-1:end])
  println("constraint2:",w_ulti[2]*betas[2,max_iters-1:end]/alphas[2,max_iters-1:end])
end

opt_alpha = zeros(No_MNO)
opt_beta = zeros(No_MNO)
opt_gamma= zeros(No_MNO)

function centralized_solver()
  prob = Model(solver=IpoptSolver(tol=1e-7, max_iter=10000, print_level =0))
    # prob = Model(solver=IpoptSolver())
  # prob = Model(solver=BonminNLSolver())
  # prob = Model(solver=OsilBonminSolver())

  @variable(prob, beta[1:No_MNO] >= 0)
  @variable(prob, alpha[1:No_MNO]>= 0)
  @variable(prob, gamma[1:No_MNO]>= 1e-8 )
  @variable(prob, lamb_miss[1:No_MNO]>=0)

  @NLobjective(prob, Min,  sum(w_ulti[i]*beta[i]/alpha[i] for i=1:No_MNO) + w_tradeoff*( sum(lamb_miss[i]/(gamma[i]*B-lamb_miss[i]) for i=1:No_MNO) ) )

  @constraint(prob, sum(beta[i] for i =1:No_MNO) == 1 )
  @constraint(prob, sum(alpha[i] for i =1:No_MNO) == 1 )
  @constraint(prob, sum(gamma[i] for i =1:No_MNO) == 1 )

  for i = 1:No_MNO
    @constraint(prob,  w_ulti[i]*beta[i]  - (1 - epsilon)*alpha[i] <= 0 )
    @constraint(prob, lamb_miss[i] +1e-8 <= gamma[i]*B )
    @constraint(prob, lamb_miss[i] == lamb_n[i] - w_lamb_hit_n[i] * beta[i])
  end


  status = solve(prob)
  println("Solve Status: ",status)
  opt_alpha[:] = getvalue(alpha)[:]
  opt_beta[:] = getvalue(beta)[:]
  opt_gamma[:] = getvalue(gamma)[:]
  global opt_cost = cal_cost_MVNOs(opt_alpha,opt_beta,opt_gamma)

  println("alpha:",opt_alpha)
  println("beta:" ,opt_beta)
  println("gamma:",opt_gamma)
  println("total_cost",opt_cost);
end

centralized_solver()

if(FLOT_ALL_FIGS)
  read_results()
  plt_all_figures()
else
  # interval=0.005
  # exhaustive_search(interval)
  main()
end
