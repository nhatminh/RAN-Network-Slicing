fig_size = (7.1,4.8)
label_fontsize = 17
patterns = ["","."]

function plt_figures()
  colors=["b","r","k"]
  figure(1,figsize=fig_size)
  plot(total_cost,color=colors[1],label="Scheme $ALG_MODE")
  plot(opt_cost*ones(max_iters),color=colors[2],linestyle="--",label="Exhaustive Search")
  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("Total Cost",fontsize=label_fontsize)
  legend(loc="best")
  title("Total Cost Convergence")
  savefig( string("figs2/cost_Convergence_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(2,figsize=fig_size)
  for i=1:No_MNO
    plot(alphas[i,:],color=colors[1],linestyle="-",label="alpha")
    plot(betas[i,:] ,color=colors[1],linestyle="-",label="beta")
    plot(gammas[i,:] ,color=colors[1],linestyle="-",label="gamma")
  end
  for i=1:No_MNO
    plot(opt_alpha[i]*ones(max_iters),color=colors[2],linestyle="--",label="alpha")
    plot(opt_beta[i] *ones(max_iters),color=colors[2],linestyle="--",label="beta")
    plot(opt_gamma[i]*ones(max_iters) ,color=colors[2],linestyle="--",label="gamma")
  end
  # plot(alphas[1,:],color=colors[1],linestyle="--",label="alpha1")
  # plot(alphas[2,:],color=colors[1],linestyle="--",label="alpha2")
  # plot(betas[1,:],color=colors[2],linestyle="--",label="beta1")
  # plot(betas[2,:],color=colors[2],linestyle="--",label="beta2")
  # plot(gammas[1,:],color=colors[3],linestyle="--",label="gamma1")
  # plot(gammas[2,:],color=colors[3],linestyle="--",label="gamma2")
  # legend(loc="best")
  title("alpha,beta,gamma")
  savefig( string("figs2/alpha_beta_gamma_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(3,figsize=fig_size)
  plot(s_fes[1,:],color=colors[1])
  plot(s_fes[2,:],color=colors[2])
  legend(loc="best")
  title("s feasible")
  savefig( string("figs2/slack_variables_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(4,figsize=fig_size)
  plot(dual_vars1,color=colors[1],label="dual_beta")
  plot(dual_vars2,color=colors[2],label="dual_alpha")
  plot(dual_vars3,color=colors[3],label="dual_gamma")
  legend(loc="best")
  title("dual")
  savefig( string("figs2/dual_variables_",ALG_MODE,"_",No_MNO,".pdf") )

  # plot(langragians,color=colors[2],linestyle="--",label="langragian")
end

function plt_all_figures()
  marker_size=6
  stride = 5

  colors=["m","b","k","g","r","coral"]
  algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
  markers = ["s","o","d"]
  figure(1,figsize=fig_size)
  plot(opt_cost1*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--",label=algs[6])
  plot(opt_cost*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(total_cost1[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",label=algs[1])
  # plot(total_cost2[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--",label=algs[2])
  plot(total_cost3[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="-",label=algs[3])
  # plot(total_cost4[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--",label=algs[4])


  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("Total Cost",fontsize=label_fontsize)
  legend(loc="best")
  # title("Total Cost Convergence")
  tight_layout()
  savefig( string("figs2/cost_Convergence_all_",No_MNO,".pdf") )

  figure(2,figsize=fig_size)
  plot(opt_alpha1[1]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--",label=algs[6])
  plot(opt_alpha1[2]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="-")
  plot(opt_alpha[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(opt_alpha[2]*ones(max_iters),color=colors[5],linestyle=":")
  # plot(alphas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # plot(alphas1[2,:],color=colors[1],linestyle="--")
  plot(alphas2[1,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--",label=algs[2])
  plot(alphas2[2,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--")
  # plot(alphas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # plot(alphas3[2,:],color=colors[3],linestyle="--")
  plot(alphas4[1,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--",label=algs[4])
  plot(alphas4[2,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--")



  legend(loc="best")
  # title("alpha")
  tight_layout()
  savefig( string("figs2/alpha_all_",No_MNO,".pdf") )

  figure(3,figsize=fig_size)
  plot(opt_beta1[1]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--",label=algs[6])
  plot(opt_beta1[2]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--")
  plot(opt_beta[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(opt_beta[2]*ones(max_iters),color=colors[5],linestyle=":")
  # plot(betas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # plot(betas1[2,:],color=colors[1],linestyle="--")
  plot(betas2[1,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--",label=algs[2])
  plot(betas2[2,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--")
  # plot(betas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # plot(betas3[2,:],color=colors[3],linestyle="--")
  plot(betas4[1,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--",label=algs[4])
  plot(betas4[2,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--")


  legend(loc="best")
  # title("beta")
  tight_layout()
  savefig( string("figs2/beta_all_",No_MNO,".pdf") )

  figure(4,figsize=fig_size)
  plot(opt_gamma1[1]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--",label=algs[6])
  plot(opt_gamma1[2]*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--")
  plot(opt_gamma[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(opt_gamma[2]*ones(max_iters),color=colors[5],linestyle=":")

  # plot(gammas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # plot(gammas1[2,:],color=colors[1],linestyle="--")
  plot(gammas2[1,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--",label=algs[2])
  plot(gammas2[2,1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--")
  # plot(gammas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # plot(gammas3[2,:],color=colors[3],linestyle="--")
  plot(gammas4[1,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--",label=algs[4])
  plot(gammas4[2,1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--")

  legend(loc="best")
  # title("gamma")
  tight_layout()
  savefig( string("figs2/gamma_all_",No_MNO,".pdf") )
end

function read_results()
  h5open("figs2/results1.h5", "r") do file
     global alphas1 = read(file, "alphas")
     global betas1 = read(file, "betas")
     global gammas1 = read(file, "gammas")
     global total_cost1 = read(file, "total_cost")
     global opt_alpha1 = read(file, "opt_alpha")
     global opt_beta1  = read(file, "opt_beta")
     global opt_gamma1 = read(file, "opt_gamma")
     global opt_cost1  = read(file, "opt_cost")

  end
  h5open("figs2/results2.h5", "r") do file
     global alphas2 = read(file, "alphas")
     global betas2  = read(file, "betas")
     global gammas2 = read(file, "gammas")
     global total_cost2 = read(file, "total_cost")
  end
  h5open("figs2/results3.h5", "r") do file
     global alphas3 = read(file, "alphas")
     global betas3  = read(file, "betas")
     global gammas3 = read(file, "gammas")
     global total_cost3 = read(file, "total_cost")
  end
  h5open("figs2/results4.h5", "r") do file
     global alphas4 = read(file, "alphas")
     global betas4  = read(file, "betas")
     global gammas4 = read(file, "gammas")
     global total_cost4 = read(file, "total_cost")
  end

end
