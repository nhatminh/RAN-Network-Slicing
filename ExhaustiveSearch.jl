opt_alpha = zeros(No_MNO)
opt_beta  = zeros(No_MNO)
opt_gamma = zeros(No_MNO)

opt_alpha1 = [0.67475,0.32525]
opt_beta1 =  [0.621,0.379]
opt_gamma1 = [0.0,1.0]
opt_cost1 = 2.198435552190618

function exhaustive_search(interval)
  s_alphas = zeros(No_MNO)
  s_betas  = zeros(No_MNO)
  s_gammas = zeros(No_MNO)
  global opt_cost1=1000

  for beta1 =0: interval:1
    s_betas[2] = 1-beta1
    for alpha1 = (w_ulti[1]*beta1/(1-epsilon)): interval:1
      s_alphas[2] = 1-alpha1
      for gamma1 = 0: interval:1
          lamb_miss1 = max(0,lamb_n[1] - w_lamb_hit_n[1]*beta1)
          if (lamb_miss1 < gamma1*B || (gamma1 == 0 && lamb_miss1 == 0))  #lambda_miss < gamma*B
            s_gammas[2] = 1- gamma1
            if (check_feasibility(s_alphas[2],s_betas[2],s_gammas[2], 2))
              s_alphas[1] = alpha1
              s_betas[1]  = beta1
              s_gammas[1] = gamma1
              total_cost = cal_cost_MVNOs(s_alphas, s_betas, s_gammas)
              if (total_cost < opt_cost1)
                opt_alpha1[:] = s_alphas[:]
                opt_beta1[:]  = s_betas[:]
                opt_gamma1[:] = s_gammas[:]
                opt_cost1 = total_cost
              end
            end
          end
      end
    end
  end
  println("alpha: ",opt_alpha1)
  println("beta: " ,opt_beta1)
  println("gamma: ",opt_gamma1)
  println("Opt_cost: ",opt_cost1)
end

function check_feasibility(alpha, beta, gamma, op_idx)
  lamb_miss = max(0,lamb_n[op_idx] - w_lamb_hit_n[op_idx]*beta)
  fes_Flag = true
  if (alpha == 0) || (alpha < w_ulti[2]*beta/(1-epsilon))
    fes_Flag = false
  elseif (lamb_miss >= gamma*B) && (gamma > 0 || lamb_miss > 0)
      fes_Flag= false
  end
  return fes_Flag
end

function check_feasibility1(alphas, betas, gammas)
  lamb_miss = max(0,lamb_n - w_lamb_hit_n.*betas)
  fes_Flag = true
  for i = 1:n
    if (alphas[i] == 0) || (alphas[i] < w_ulti[i]*betas[i]/(1-epsilon))
      fes_Flag = false
      break
    elseif (lamb_miss[i] > gammas[i]*B) && (gammas[i] > 0 || lamb_miss[i] > 0)
      fes_Flag= false
      break
    end
  end
  return fes_Flag
end

function cal_cost_MVNOs(s_alphas, s_betas, s_gammas)
  lamb_miss = max(0,lamb_n - w_lamb_hit_n .*s_betas)
  backhaul_cost = 0

  for i = 1:No_MNO
    if (lamb_miss[i] < s_gammas[i]*B ) #zero cost when gamma*B = lambda_miss = 0
      backhaul_cost += lamb_miss[i]/(s_gammas[i]*B - lamb_miss[i])
    end
  end

  total_cost =  vecdot(w_ulti,s_betas./s_alphas) + w_tradeoff*backhaul_cost
  return total_cost
end
