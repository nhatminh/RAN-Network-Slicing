# Let us first make the Convex and MultiConvex modules available
using Convex
using ECOS
using Gurobi
using SCS
using Ipopt
using PyPlot

n=2
max_iters = 400
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
w_init   = 1e-1
w_ulti = [0.9, 0.8 ]
rho = 0.5
rho1 = 10.
w_tradeoff =1

function centralized_solver(op::MOperator, alpha = 2)
  prob = Model(solver=IpoptSolver(tol=5e-6, max_iter=10000, print_level =2))
  # prob = Model(solver=BonminNLSolver())
  # prob = Model(solver=OsilBonminSolver())

  @variable(prob, beta[1:No_MNO] >= 0)
  @variable(prob, alpha[1:No_MNO]>0)
  @variable(prob, gamma[1:No_MNO]>=0)
  @variable(prob, lambda_miss[1:Numb_BS]>=0)

  @NLobjective(prob, Min,  sum{w_ulti[i]*beta[i]/alpha[i], i=1:No_MNO} + w_tradeoff*( sum{max(0,gamma[i]*B/(gamma[i]*B-lamb_miss[i])),i=1:No_MNO} ) )

  @constraint(prob, sum{beta[i], i =1:No_MNO} == 1 )
  @constraint(prob, sum{alpha[i], i =1:No_MNO} == 1 )
  @constraint(prob, sum{gamma[i], i =1:No_MNO} == 1 )
  for i = 1:No_MNO
    @constraint(prob,  w_ulti[i]*beta  - (1 - epsilon)*alpha[i] <= 0 )
    @constraint(prob, lambda_miss[i] < gamma[i]*B )
  end


  status = solve(prob)
  println("Dis Solve: ",status)
  println(getvalue(alpha))
end

function beta_update0(k)
  beta = Variable(n)

  obj  =  vecdot(w_ulti,beta./alphas[:,k]) + rho*sumsquares(beta - betas[:,k])+
          # sum(beta)
          sum( max( beta,3*beta -2/3*gammas[:,k],10*beta - 16/3*gammas[:,k],70*beta - 178/3*gammas[:,k],
          500*beta - 1468/3*gammas[:,k],5000*beta - 16318/3*gammas[:,k] ) )
  problem = minimize(obj)
  problem.constraints +=  [sum(beta) == 1, beta>=0, w_ulti.*beta  - (1 - epsilon)*alphas[:,k] <= 0]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)
  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  return beta.value
end

function alpha_gamma_update0(k)
  alpha = Variable(n)
  gamma = Variable(n)

  obj  = vecdot(w_ulti,betas[:,k+1].*invpos(alpha)) + rho*sumsquares(alpha - alphas[:,k]) + rho*sumsquares(gamma - gammas[:,k]) +
  # sum(-gamma)
  sum( max( betas[:,k],3*betas[:,k] -2/3*gamma,10*betas[:,k] - 16/3*gamma,70*betas[:,k] - 178/3*gamma,
  500*betas[:,k] - 1468/3*gamma,5000*betas[:,k] - 16318/3*gamma ) )
  problem = minimize(obj)
  problem.constraints +=  [sum(alpha) == 1,sum(gamma) == 1, alpha >= 0, gamma >= 0, w_ulti.*betas[:,k+1] - (1 - epsilon)*alpha <=0]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("alphagamma-subprob: ",problem.status)
  # println("alpha:", alpha.value)
  # println("gamma:", gamma.value)
  return alpha.value, gamma.value, [0,0]
end

function beta_i_update0(k,i)
  beta = Variable()
  sub_obj = beta
  for idx = 1:n
      if idx< i
        sub_obj += betas[idx,k+1]
      elseif idx > i
        sub_obj += betas[idx,k]
      end
  end

  obj  =  w_ulti[i]*beta/alphas[i,k] + rho*square(beta - betas[i,k]) + rho1*square(sub_obj - 1 - dual_vars1[k]) +
          # beta
          max( beta,3*beta -2/3*gammas[i,k],10*beta - 16/3*gammas[i,k],70*beta - 178/3*gammas[i,k],
          500*beta - 1468/3*gammas[i,k],5000*beta - 16318/3*gammas[i,k] )

  problem = minimize(obj)
  problem.constraints +=  [beta>=0, w_ulti[i]*beta  - (1 - epsilon)*alphas[i,k] <=0]

  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_i_update0(k,i)
  alpha = Variable()
  gamma = Variable()

  sub_obj2 = alpha
  sub_obj3 = gamma
  for idx = 1:n
      if idx< i
        sub_obj2 += alphas[idx,k+1]
        sub_obj3 += gammas[idx,k+1]
      elseif idx > i
        sub_obj2 += alphas[idx,k]
        sub_obj3 += gammas[idx,k]
      end
  end

  obj  = w_ulti[i]*betas[i,k+1]/(alpha)+ rho*((alpha - alphas[i,k])^2)+ rho*((gamma - gammas[i,k])^2) +
  rho1*((sub_obj2 - 1 - dual_vars2[k])^2) + rho1*((sub_obj3 - 1 - dual_vars3[k])^2) +
  # (-gamma)
  max( betas[i,k],3*betas[i,k] -2/3*gamma,10*betas[i,k] - 16/3*gamma,70*betas[i,k] - 178/3*gamma,
  500*betas[i,k] - 1468/3*gamma,5000*betas[i,k] - 16318/3*gamma )

  problem = minimize(obj)
  problem.constraints +=  [alpha >= 0, gamma >= 0, w_ulti[i]*betas[i,k+1] - (1 - epsilon)*alpha <=0]
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

function beta_i_update(k,i)
  beta = Variable()
  s = Variable()
  sub_obj = beta
  for idx = 1:n
      if idx< i
        sub_obj += betas[idx,k+1]
      elseif idx > i
        sub_obj += betas[idx,k]
      end
  end

  obj  =  w_ulti[i]*beta/alphas[i,k] +w*s+ rho*square(beta - betas[i,k]) + rho1*square(sub_obj - 1 - dual_vars1[k]) +
          # beta
          max( beta,3*beta -2/3*gammas[i,k],10*beta - 16/3*gammas[i,k],70*beta - 178/3*gammas[i,k],
          500*beta - 1468/3*gammas[i,k],5000*beta - 16318/3*gammas[i,k] )

  problem = minimize(obj)
  problem.constraints +=  [beta>=0, w_ulti[i]*beta  - (1 - epsilon)*alphas[i,k] <=s, s>=0]

  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_i_update(k,i)
  alpha = Variable()
  gamma = Variable()
  s = Variable()

  sub_obj2 = alpha
  sub_obj3 = gamma
  for idx = 1:n
      if idx< i
        sub_obj2 += alphas[idx,k+1]
        sub_obj3 += gammas[idx,k+1]
      elseif idx > i
        sub_obj2 += alphas[idx,k]
        sub_obj3 += gammas[idx,k]
      end
  end

  obj  = w_ulti[i]*betas[i,k+1]/(alpha)+ w*s + rho*((alpha - alphas[i,k])^2)+ rho*((gamma - gammas[i,k])^2) +
  rho1*((sub_obj2 - 1 - dual_vars2[k])^2) + rho1*((sub_obj3 - 1 - dual_vars3[k])^2) +
  # (-gamma)
  max( betas[i,k],3*betas[i,k] -2/3*gamma,10*betas[i,k] - 16/3*gamma,70*betas[i,k] - 178/3*gamma,
  500*betas[i,k] - 1468/3*gamma,5000*betas[i,k] - 16318/3*gamma )

  problem = minimize(obj)
  problem.constraints +=  [alpha >= 0, gamma >= 0, w_ulti[i]*betas[i,k+1] - (1 - epsilon)*alpha <=s, s>=0]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  return alpha.value, gamma.value, s.value
end

function beta_update(k)
  beta = Variable(n)
  s = Variable(n)

  obj  =  vecdot(w_ulti,beta./alphas[:,k]) + w*sum(s) + rho*sumsquares(beta - betas[:,k]) +
          # sum(beta)
          sum( max( beta,3*beta -2/3*gammas[:,k],10*beta - 16/3*gammas[:,k],70*beta - 178/3*gammas[:,k],
          500*beta - 1468/3*gammas[:,k],5000*beta - 16318/3*gammas[:,k] ) )

  problem = minimize(obj)
  problem.constraints +=  [sum(beta) == 1, beta>=0, w_ulti.*beta  - (1 - epsilon)*alphas[:,k] <= s, s>= 0]
  # solve!(prob, GurobiSolver(),verbose=false)
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)
  # solve!(prob, SCSSolver(verbose=false),verbose=false)
  # println("beta-subprob: ",problem.status)
  # println("beta:", beta.value)
  # println("s:", s.value)
  return beta.value
end

function alpha_gamma_update(k)
  alpha = Variable(n)
  gamma = Variable(n)
  s = Variable(n)

  obj  = vecdot(w_ulti,betas[:,k+1].*invpos(alpha)) + w*sum(s) + rho*sumsquares(alpha - alphas[:,k]) + rho*sumsquares(gamma - gammas[:,k]) +
  # sum(-gamma)
  sum( max( betas[:,k],3*betas[:,k] -2/3*gamma,10*betas[:,k] - 16/3*gamma,70*betas[:,k] - 178/3*gamma,
  500*betas[:,k] - 1468/3*gamma,5000*betas[:,k] - 16318/3*gamma ) )

  problem = minimize(obj)
  problem.constraints +=  [sum(alpha) == 1,sum(gamma) == 1, alpha >= 0, gamma >= 0, w_ulti.*betas[:,k+1] - (1 - epsilon)*alpha <= s, s>= 0]
  solve!(problem, ECOSSolver(verbose=false, max_iters = 1000),verbose=false)

  # println("alphagamma-subprob: ",problem.status)
  # println("alpha:", alpha.value)
  # println("gamma:", gamma.value)
  # println("s:", s.value)
  return alpha.value, gamma.value, s.value
end

function cal_total_cost(k)
  total_cost[k] = sum(betas[:,k]/alphas[:,k]) + sum(betas[:,k]) - sum(gammas[:,k])
end

function cal_langragian(k)
  langragians[k] = sum(betas[:,k]/alphas[:,k]) + sum(betas[:,k]) - sum(gammas[:,k]) + rho*(norm(alphas[:,k] - alphas[:,k-1])^2)
   + rho*(norm(gammas[:,k] - gammas[:,k-1])^2) + rho*(norm(betas[:,k] - betas[:,k-1])^2)
end

function plt_figures()
  colors=["b","r","k"]
  figure(1)
  plot(total_cost,color=colors[1],label="total cost")
  legend(loc="best")
  title("Total Cost")

  figure(2)
  plot(alphas[1,:],color=colors[1],linestyle="--",label="alpha1")
  plot(alphas[2,:],color=colors[1],linestyle="--",label="alpha2")
  plot(betas[1,:],color=colors[2],linestyle="--",label="beta1")
  plot(betas[2,:],color=colors[2],linestyle="--",label="beta2")
  plot(gammas[1,:],color=colors[3],linestyle="--",label="gamma1")
  plot(gammas[2,:],color=colors[3],linestyle="--",label="gamma2")
  legend(loc="best")
  title("alpha,beta,gamma")
  figure(3)
  plot(s_fes[1,:],color=colors[1])
  plot(s_fes[2,:],color=colors[2])
  legend(loc="best")
  title("s feasible")

  figure(4)
  plot(dual_vars1,color=colors[1],label="dual_beta")
  plot(dual_vars2,color=colors[2],label="dual_alpha")
  plot(dual_vars3,color=colors[3],label="dual_gamma")
  legend(loc="best")
  title("dual")

  # plot(langragians,color=colors[2],linestyle="--",label="langragian")
end

function main()
  alphas[:,1] = [0.5 0.5]
  betas[:,1]  = [0.5 0.5]
  gammas[:,1] = [0.5 0.5]
  cal_total_cost(1)
  global rho = 0.5
  global w = w_init
  # cal_langragian(1)

  for k=1:max_iters-1
    rho = rho
    w = min(max_w,w *w_factor)
    # betas[:,k+1]  = beta_update0(k)
    # alphas[:,k+1], gammas[:,k+1], s_fes[:,k+1] = alpha_gamma_update0(k)
    betas[:,k+1]  = beta_update(k)
    alphas[:,k+1], gammas[:,k+1], s_fes[:,k+1] = alpha_gamma_update(k)
    cal_total_cost(k+1)
    # cal_langragian(k+1)
  end
  plt_figures()

  # for k=1:max_iters-1
  #   # println("---- Iteration: ",k)
  #   rho = rho
  #   w = min(max_w,w *w_factor)
  #   for i =1:n
  #     betas[i,k+1]  = beta_i_update0(k,i)
  #     alphas[i,k+1], gammas[i,k+1], s_fes[i,k+1] = alpha_gamma_i_update0(k,i)
  #     dual_vars1[k+1] = dual_vars1[k] - sum(betas[:,k+1])  + 1
  #     dual_vars2[k+1] = dual_vars2[k] - sum(alphas[:,k+1]) + 1
  #     dual_vars3[k+1] = dual_vars3[k] - sum(gammas[:,k+1]) + 1
  #   end
  #   cal_total_cost(k+1)
  #   # cal_langragian(k+1)
  # end
  # plt_figures()

  println("alpha",alphas[:,max_iters-2:max_iters])
  println("beta",betas[:,max_iters-2:max_iters])
  println("gamma",gammas[:,max_iters-2:max_iters])
  println("total_cost",total_cost[max_iters-2:max_iters])
  println("s_fes:",s_fes[:,max_iters-2:max_iters])
  println("constraint1:",w_ulti*betas[1,max_iters]/alphas[1,max_iters])
  println("constraint2:",w_ulti*betas[2,max_iters]/alphas[2,max_iters])
end

main()
