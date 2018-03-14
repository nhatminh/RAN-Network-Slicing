using Convex
using ECOS
using Gurobi
using SCS
using PyPlot
using HDF5
using JuMP
using Ipopt

REUSED_INPUT = true
FLOT_ALL_FIGS = false
F= 5*10^4
S= 220
B= 500
L= 100*1024*8  #File size
No_MNO=2
max_iters = 55
ALG_MODE = 3 #1, 2, 3, 4
w_tradeoff= 1.
# F=5*10^4
# S=100
# B=400
# L= 150*1024*8  #File size
# No_MNO=2
# max_iters = 3500
# ALG_MODE = 3 #1, 2, 3, 4
# w_tradeoff= 1e-2

input_filename = "weights2.h5"
if (REUSED_INPUT)
   h5open(input_filename, "r") do file
      global lamb_n = read(file, "lamb_n")
      global w_ulti = read(file, "w_ulti")
      global w_lamb_hit_n = read(file, "w_lamb_hit_n")
      global N_UE_s = read(file, "N_UE_s")
   end
else
  include("Setting$No_MNO-MC.jl")
  random_input()
  # include("Setting$No_MNO.jl")
  h5open(input_filename, "w") do file
      write(file, "lamb_n", lamb_n)
      write(file, "w_ulti", w_ulti)
      write(file, "w_lamb_hit_n", w_lamb_hit_n)
      write(file, "N_UE_s",N_UE_s)
  end
end

include("ExhaustiveSearch.jl")
include("Plot_Fig2.jl")

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

# rho = 0.05 * 2
# rho1 = 10. * 2
rho = 8
rho1 = 12
Jacobian_step = 1.
if ALG_MODE > 2
  delta = 1.
  rho = (rho1*(No_MNO/(2-Jacobian_step)-1) + delta)
end

function beta_obj(beta,lamb_miss, k)
  obj  =  vecdot(w_ulti,beta./alphas[:,k]) + rho/2*sumsquares(beta - betas[:,k]) +
  # w_tradeoff*(-No_MNO + sum((gammas[:,k]*B).*invpos(gammas[:,k]*B-lamb_miss)) )
  w_tradeoff*(-No_MNO + sum( max(0,(gammas[:,k]*B).*invpos(gammas[:,k]*B-lamb_miss)))  )

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

  return obj
end


function alpha_gamma_obj(alpha,gamma,lamb_miss, k)
  obj  = vecdot(w_ulti,betas[:,k+1].*invpos(alpha)) + rho/2*sumsquares(alpha - alphas[:,k]) + rho/2*sumsquares(gamma - gammas[:,k]) +
  # w_tradeoff*sum(lamb_miss.*invpos(gamma*B-lamb_miss))
  w_tradeoff*sum( max(0,lamb_miss.*invpos(gamma*B-lamb_miss)) )

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
  end
end

function save_results()
  h5open("figs2/results$ALG_MODE.h5", "w") do file
     write(file, "alphas", alphas)
     write(file, "betas" , betas)
     write(file, "gammas", gammas)
     write(file, "total_cost", total_cost)
     write(file, "opt_alpha" , opt_alpha1)
     write(file, "opt_beta"  , opt_beta1)
     write(file, "opt_gamma" , opt_gamma1)
     write(file, "opt_cost"  , opt_cost1)
     write(file, "opt_solver"  , opt_cost)
  end
end

function main()
  # alphas[:,1] = 1/No_MNO *ones(No_MNO)
  # betas[:,1]  = 1/No_MNO *ones(No_MNO)
  # gammas[:,1] = 1/No_MNO *ones(No_MNO)
  alphas[:,1] = 0.3 *ones(No_MNO)
  betas[:,1]  = 0.3 *ones(No_MNO)
  gammas[:,1] = 0.3 *ones(No_MNO)

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
  interval=0.005
  exhaustive_search(interval)
  main()
end
