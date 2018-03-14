using Convex
using ECOS
using Gurobi
using SCS
using Ipopt
using PyPlot
using HDF5
using JuMP

FLOT_MC_FIGS = true
MONTE_CARLO = true
REUSED_INPUT = false

F= 5*10^4
S= 220
B= 1500
L= 100*1024*8  #File size
No_MNO=4
max_iters = 200
ALG_MODE = 1 #1, 2, 3, 4
w_tradeoff= 1.
N_SIMs = 300
Cost_MC=zeros(4,3,N_SIMs)
Converged_Iters_MC=zeros(4,N_SIMs)
Cost_MC_Cen1=zeros(3,N_SIMs)
Cost_MC_Cen2=zeros(3,N_SIMs)
Cost_MC_Cen3=zeros(3,N_SIMs)
Cost_MC_Cen3_e1=zeros(2,No_MNO,N_SIMs)
Cost_MC_Cen3_e2=zeros(2,No_MNO,N_SIMs)

include("Setting$No_MNO-MC.jl")
include("Plot_MC.jl")
include("Plot_Fig_general.jl")

n=No_MNO
epsilon = 1e-3
STOP_EPS =1e-2

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
                          lamb_miss + 1e-8*ones(n) <= gammas[:,k]*B]
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
  problem.constraints +=  [sum(alpha) == 1,sum(gamma) == 1,alpha >= 1e-8,gamma >= 0, w_ulti.*betas[:,k+1] - (1 - epsilon)*alpha <=0,
                          lamb_miss + 1e-8*ones(n) <= gamma*B]
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

function cal_total_cost(k, s)
  lamb_miss = max(0,lamb_n - w_lamb_hit_n .*betas[:,k])
  backhaul_cost = sum(lamb_miss./(gammas[:,k]*B - lamb_miss) )
  backhaul_cost = 0

  for i=1:No_MNO
    if gammas[i,k]*B <= lamb_miss[i]
      println("M/M/1 Violation")
    else
      backhaul_cost += lamb_miss[i]/(gammas[i,k]*B - lamb_miss[i])
    end
  end

  total_cost[k] =  vecdot(w_ulti,betas[:,k]./alphas[:,k]) + w_tradeoff*backhaul_cost
  ## LOG for each simulation result
  # Cost_MC[ALG_MODE,1,s] = vecdot(w_ulti,betas[:,k]./alphas[:,k])
  # Cost_MC[ALG_MODE,2,s] = backhaul_cost
  # Cost_MC[ALG_MODE,3,s] = total_cost[k]
  Converged_Iters_MC[ALG_MODE,s] = k+1

  if(k==max_iters)
    println("Test lamb_miss:", lamb_miss)
    println("Test utilization:",vecdot(w_ulti,betas[:,k]./alphas[:,k]))
    println("Test miss cost:", lamb_miss./(gammas[:,k]*B - lamb_miss) )
  end
end

function save_results()
  h5open(string(folder,"results$ALG_MODE.h5"), "w") do file
     write(file, "alphas", alphas)
     write(file, "betas" , betas)
     write(file, "gammas", gammas)
     write(file, "total_cost", total_cost)
  end
end

function save_MC_results()
  h5open(string(folder,"MC$N_SIMs.h5"), "w") do file
     write(file, "Cost_MC", Cost_MC)
     write(file, "Converged_Iters_MC" , Converged_Iters_MC)
  end
end

function save_MC_Cen_results()
  h5open(string(folder,"MC$N_SIMs-Cen1.h5"), "w") do file
     write(file, "Cost_MC_Cen1", Cost_MC_Cen1)
     write(file, "Cost_MC_Cen2", Cost_MC_Cen2)
     write(file, "Cost_MC_Cen3", Cost_MC_Cen3)
     write(file, "Converged_Iters_MC" , Converged_Iters_MC)
  end
end

function save_MC_Cen_results_eps()
  h5open(string(folder,"MC$N_SIMs-Cen-eps.h5"), "w") do file
     write(file, "Cost_MC_Cen3_e1", Cost_MC_Cen3_e1) #epsilon = 0.001
     write(file, "Cost_MC_Cen3_e2", Cost_MC_Cen3_e2)  #epsilon = 0.1
  end
end

function main(s)
  # alphas[:,1] = 1/No_MNO *ones(No_MNO)
  # betas[:,1]  = 1/No_MNO *ones(No_MNO)
  # gammas[:,1] = 1/No_MNO *ones(No_MNO)
  alphas[:,1] = 0.1*ones(No_MNO)
  betas[:,1]  = 0.1*ones(No_MNO)
  gammas[:,1] = 0.1*ones(No_MNO)
  cal_total_cost(1,s)
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
    cal_total_cost(k+1,s)
    if abs(total_cost[k+1]-total_cost[k]) < STOP_EPS
      break;
    end
  end

  if(Converged_Iters_MC[ALG_MODE,s] == max_iters+1)
    h5open("Weights4_NonConverged", "w") do file
        write(file, "lamb_n", lamb_n)
        write(file, "w_ulti", w_ulti)
        write(file, "w_lamb_hit_n", w_lamb_hit_n)
    end
  end

  if(REUSED_INPUT)
    plt_figures()
  end

  # println("alpha",alphas[:,end])
  # println("alpha_opt",opt_alpha)
  # println("beta",betas[:,end])
  # println("beta_opt",opt_beta)
  # println("gamma",gammas[:,end])
  # println("gamma_opt",opt_gamma)
  # println("total_cost:",total_cost[end])
end

opt_alpha = zeros(No_MNO)
opt_beta = zeros(No_MNO)
opt_gamma= zeros(No_MNO)

function cal_cost_MVNOs(s_alphas, s_betas, s_gammas, s)
  lamb_miss = max(0,lamb_n - w_lamb_hit_n .*s_betas)
  backhaul_cost = 0

  for i = 1:No_MNO
    if (lamb_miss[i] < s_gammas[i]*B ) #zero cost when gamma*B = lambda_miss = 0
      backhaul_cost += lamb_miss[i]/(s_gammas[i]*B - lamb_miss[i])
      # println(lamb_n[i] - w_lamb_hit_n[i] *s_betas[i])
    end
  end

  total_cost =  vecdot(w_ulti,s_betas./s_alphas) + w_tradeoff*backhaul_cost
  # if(w_tradeoff == 0.01)
  #   Cost_MC_Cen1[1,s] = vecdot(w_ulti,s_betas./s_alphas)
  #   Cost_MC_Cen1[2,s] = backhaul_cost
  #   Cost_MC_Cen1[3,s] = total_cost
  # elseif (w_tradeoff == 0.1)
  #   Cost_MC_Cen2[1,s] = vecdot(w_ulti,s_betas./s_alphas)
  #   Cost_MC_Cen2[2,s] = backhaul_cost
  #   Cost_MC_Cen2[3,s] = total_cost
  # else
  #   Cost_MC_Cen3[1,s] = vecdot(w_ulti,s_betas./s_alphas)
  #   Cost_MC_Cen3[2,s] = backhaul_cost
  #   Cost_MC_Cen3[3,s] = total_cost
  # end
  if(epsilon == 0.001)
    for i = 1:No_MNO
      Cost_MC_Cen3_e1[1,i,s] = w_ulti[i] * s_betas[i]/s_alphas[i]
      if (lamb_miss[i] < s_gammas[i]*B )
        Cost_MC_Cen3_e1[2,i,s] = lamb_miss[i]/(s_gammas[i]*B - lamb_miss[i])
      end
    end
    # println("Cost_MC_Cen3_e1:",Cost_MC_Cen3_e1[1,:,s])
    # println("Cost_MC_Cen3_e1:",Cost_MC_Cen3_e1[2,:,s])
  elseif epsilon == 0.1
    for i = 1:No_MNO
      Cost_MC_Cen3_e2[1,i,s] = w_ulti[i] * s_betas[i]/s_alphas[i]
      if (lamb_miss[i] < s_gammas[i]*B )
        Cost_MC_Cen3_e2[2,i,s] = lamb_miss[i]/(s_gammas[i]*B - lamb_miss[i])
      end
    end
    # println("Cost_MC_Cen3_e2:",Cost_MC_Cen3_e2[1,:,s])
    # println("Cost_MC_Cen3_e2:",Cost_MC_Cen3_e2[2,:,s])
  end

  # println("Utilization2:",vecdot(w_ulti,s_betas./s_alphas))
  # println("backhaul_cost2:",backhaul_cost)
  # println("totalcost:",total_cost)
  return total_cost
end

function centralized_solver(s)
  println(w_tradeoff)
  prob = Model(solver=IpoptSolver(tol=1e-9, max_iter=10000, print_level =0))
  # prob = Model(solver=IpoptSolver())
  # prob = Model(solver=BonminNLSolver())
  # prob = Model(solver=OsilBonminSolver())

  @variable(prob, beta[1:No_MNO] >= 0)
  @variable(prob, alpha[1:No_MNO]>= 1e-10)
  @variable(prob, gamma[1:No_MNO]>= 0 )
  @variable(prob, lamb_miss[1:No_MNO]>=0)

  @NLobjective(prob, Min,  sum(w_ulti[i]*beta[i]/alpha[i] for i=1:No_MNO) +
  # w_tradeoff*( sum((lamb_n[i] - w_lamb_hit_n[i] * beta[i])/(gamma[i]*B-(lamb_n[i] - w_lamb_hit_n[i] * beta[i])) for i=1:No_MNO) ) )
  w_tradeoff*( sum(lamb_miss[i]/(gamma[i]*B-lamb_miss[i]) for i=1:No_MNO) ) )

  @constraint(prob, sum(beta[i] for i =1:No_MNO) == 1 )
  @constraint(prob, sum(alpha[i] for i =1:No_MNO) == 1 )
  @constraint(prob, sum(gamma[i] for i =1:No_MNO) == 1 )

  # @constraint(prob, sum(beta[i] for i =1:No_MNO) >= 1-1e-3 )
  # @constraint(prob, sum(alpha[i] for i =1:No_MNO) >= 1-1e-3 )
  # @constraint(prob, sum(gamma[i] for i =1:No_MNO) >= 1-1e-3 )

  for i = 1:No_MNO
    @constraint(prob,  w_ulti[i]*beta[i]  - (1 - epsilon)*alpha[i] <= 0 )
    @constraint(prob, (lamb_n[i] - w_lamb_hit_n[i] * beta[i]) +1e-10 <= gamma[i]*B )
    @constraint(prob, lamb_miss[i] == lamb_n[i] - w_lamb_hit_n[i] * beta[i])
    # @constraint(prob, lamb_n[i] >= w_lamb_hit_n[i] * beta[i])
  end


  status = solve(prob)
  println("Solve Status: ",status)
  opt_alpha[:] = getvalue(alpha)[:]
  opt_beta[:] = getvalue(beta)[:]
  opt_gamma[:] = getvalue(gamma)[:]
  global opt_cost = cal_cost_MVNOs(opt_alpha,opt_beta,opt_gamma,s)

  println("alpha:",opt_alpha)
  println("beta:" ,opt_beta)
  println("gamma:",opt_gamma)
  println("total_cost:",opt_cost);
end

function centralized_checking(alpha, beta, gamma)
  # @variable(prob, beta[1:No_MNO] >= 0)
  # @variable(prob, alpha[1:No_MNO]>= 1e-8)
  # @variable(prob, gamma[1:No_MNO]>= 0 )
  # @variable(prob, lamb_miss[1:No_MNO]>=0)
  lamb_miss = lamb_n - w_lamb_hit_n.*beta
  println("lamb_miss",lamb_miss)
  flag = true
  for i = 1:No_MNO
    if( w_ulti[i]*beta[i]  - (1 - epsilon)*alpha[i] > 0 )
      flag = false
      println("False 1")
    end
    if( lamb_miss[i] +1e-8 > gamma[i]*B )
      flag = false
      println("False 2")
    end
    if(  beta[i] < 0 )
      flag = false
      println("False 3")
    end
    if(  alpha[i] < 1e-8 )
      flag = false
      println("False 4")
    end
    if(gamma[i]<0 )
      flag = false
      println("False 5")
    end
    if(lamb_miss[i]<0 )
      flag = false
      println("False 6")
      println(lamb_miss[i])
    end
  end
  if(flag) println("Correct") end
end

function check_feasibility1(s_alphas, s_betas, s_gammas)
  lamb_miss = lamb_n - w_lamb_hit_n.*s_betas
  fes_Flag = true
  for i = 1:n
    if (s_alphas[i] == 0) || (s_alphas[i]*(1-epsilon) < w_ulti[i]*s_betas[i])
      fes_Flag = false
      break
    # elseif (lamb_miss[i] >= s_gammas[i]*B) && (s_gammas[i] > 0 || lamb_miss[i] > 0)
    elseif (lamb_miss[i] >= s_gammas[i]*B)
      fes_Flag= false
      break
    end
  end
  # println("sum_a:",sum(s_alphas))
  # println("sum_b:",sum(s_betas))
  # println("sum_g:",sum(s_gammas))
  centralized_checking(s_alphas,s_betas,s_gammas)
  return fes_Flag
end


if (FLOT_MC_FIGS)
  read_MC()
  plt_epsilon()
elseif (MONTE_CARLO)
    for s =1:N_SIMs
      println("----- SIM $s ----")
      #RUN RANDOM
      if(REUSED_INPUT)
        h5open("Weights4_NonConverged", "r") do file
           global lamb_n = read(file, "lamb_n")
           global w_ulti = read(file, "w_ulti")
           global w_lamb_hit_n = read(file, "w_lamb_hit_n")
        end
        centralized_solver(s)
        # centralized_solver()
        global ALG_MODE = 1
        # global ALG_MODE = 3
        main(s)
      else
        random_input()
        ##### VARY w_tradeoff
        # global w_tradeoff = 0.01
        # centralized_solver(s)
        # global w_tradeoff = 0.1
        # centralized_solver(s)
        # global w_tradeoff = 1
        # centralized_solver(s)
        #
        # global ALG_MODE = 1
        # main(s)
        # global ALG_MODE = 3
        # main(s)

        ##### VARY epsilon
        # global w_tradeoff = 1
        global w_tradeoff = 1
        global epsilon = 1e-3
        centralized_solver(s)
        global epsilon = 0.1
        centralized_solver(s)
      end
      # println("W_util:",w_ulti)
      # println("lamb_n",lamb_n)
      # println("W_lamb_hit_n:",w_lamb_hit_n)
    end
    println("----- Solution ----")
    println(string("# Simulations: $N_SIMs"))
    # println("* Iters:",Converged_Iters_MC)
    # println("Costs:",Cost_MC)
    # save_MC_results()
    # println("* Costs1:",Cost_MC_Cen1)
    # println("* Costs2:",Cost_MC_Cen2)
    # println("* Costs3:",Cost_MC_Cen3)

    save_MC_Cen_results_eps()
    plt_epsilon()
    # plt_Stats()

  # println("Feasibility_opt:",check_feasibility1(opt_alpha,opt_beta,opt_gamma))
  # println("Feasibility:",check_feasibility1(alphas[:,end],betas[:,end],gammas[:,end]))
  # println("Cost Checking:",cal_cost_MVNOs(opt_alpha,opt_beta,opt_gamma))
end
