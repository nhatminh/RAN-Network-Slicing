# fig_size = (7.,5.1)
fig_size = (6.,4.3)
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 4
patterns = ["","."]

folder = string("figs$No_MNO/")
function plt_figures()
  colors=["b","r","k"]
  figure(1,figsize=fig_size)
  plot(total_cost,color=colors[1],label="Scheme $ALG_MODE")
  plot(opt_cost*ones(max_iters),color=colors[2],linestyle="--",label="IpOpt Solver")
  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("Total Cost",fontsize=label_fontsize)
  legend(loc="best")
  title("Total Cost Convergence")
  savefig( string(folder,"cost_Convergence_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(2,figsize=fig_size)
  for i=1:No_MNO
    plot(alphas[i,:],color=colors[1],linestyle="-",label="alpha")
    plot(opt_alpha[i]*ones(max_iters) ,color=colors[2],linestyle="--",label="opt_alpha")
    plot(betas[i,:] ,color=colors[1],linestyle="-",label="beta")
    plot(opt_beta[i]*ones(max_iters) ,color=colors[2],linestyle="--",label="opt_beta")
    plot(gammas[i,:] ,color=colors[1],linestyle="-",label="gamma")
    plot(opt_gamma[i]*ones(max_iters) ,color=colors[2],linestyle="--",label="opt_gamma")
  end
  title("alpha,beta,gamma")
  savefig( string(folder,"alpha_beta_gamma_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(3,figsize=fig_size)
  for i=1:No_MNO
    plot(s_fes[i,:],color=colors[1])
  end
  legend(loc="best")
  title("s feasible")
  savefig( string(folder,"slack_variables_",ALG_MODE,"_",No_MNO,".pdf") )

  figure(4,figsize=fig_size)
  plot(dual_vars1,color=colors[1],label="dual_beta")
  plot(dual_vars2,color=colors[2],label="dual_alpha")
  plot(dual_vars3,color=colors[3],label="dual_gamma")
  legend(loc="best")
  title("dual")
  savefig( string(folder,"dual_variables_",ALG_MODE,"_",No_MNO,".pdf") )

  # plot(langragians,color=colors[2],linestyle="--",label="langragian")
end

function plt_all_figures()
  marker_size=6
  l_width=1.2
  # stride = convert(Int, max_iters/10) +1
  stride = 6

  colors=["m","b","k","g","r","coral"]
  algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
  markers = ["s","o","d"]

  figure(1,figsize=fig_size)
  plot(opt_cost*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(total_cost1[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",label=algs[1])
  # plot(total_cost2[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="--",label=algs[2])
  plot(total_cost3[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="-",label=algs[3])
  # plot(total_cost4[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="--",label=algs[4])

  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("Total Cost",fontsize=label_fontsize)
  legend(loc="best")
  # title("Total Cost Convergence")
  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  savefig( string(folder,"cost_Convergence_all_",No_MNO,".pdf") )

  figure(2,figsize=fig_size)
  plot(opt_alpha[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  for i=2:No_MNO
    plot(opt_alpha[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  # plot(alphas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # for i=2:No_MNO
  #   plot(alphas1[i,:],color=colors[1],linestyle="--")
  # end
  plot(alphas1[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[1])
  for i=2:No_MNO
    plot(alphas1[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  # plot(alphas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # for i=2:No_MNO
  #   plot(alphas3[i,:],color=colors[3],linestyle="--")
  # end
  plot(alphas3[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[3])
  for i=2:No_MNO
    plot(alphas3[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("\$\\alpha\$",fontsize=label_fontsize)
  legend(loc="best")
  # title("alpha")
  tight_layout()
  savefig( string(folder,"alpha_all_",No_MNO,".pdf") )

  figure(3,figsize=fig_size)
  plot(opt_beta[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  for i=2:No_MNO
    plot(opt_beta[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  # plot(betas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # for i=2:No_MNO
  #   plot(betas1[i,:],color=colors[1],linestyle="--")
  # end
  plot(betas1[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[1])
  for i=2:No_MNO
    plot(betas1[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  # plot(betas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # for i=2:No_MNO
  #   plot(betas3[i,:],color=colors[3],linestyle="--")
  # end
  plot(betas3[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[3])
  for i=2:No_MNO
    plot(betas3[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("\$\\beta\$",fontsize=label_fontsize)
  legend(loc="best")
  # title("beta")
  tight_layout()
  savefig( string(folder,"beta_all_",No_MNO,".pdf") )

  figure(4,figsize=fig_size)
  plot(opt_gamma[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  for i=2:No_MNO
    plot(opt_gamma[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  # plot(gammas1[1,:],color=colors[1],linestyle="--",label=algs[1])
  # for i=2:No_MNO
  #   plot(gammas1[i,:],color=colors[1],linestyle="--")
  # end
  plot(gammas1[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[1])
  for i=2:No_MNO
    plot(gammas1[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  # plot(gammas3[1,:],color=colors[3],linestyle="--",label=algs[3])
  # for i=2:No_MNO
  #   plot(gammas3[i,:],color=colors[3],linestyle="--")
  # end
  plot(gammas3[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[3])
  for i=2:No_MNO
    plot(gammas3[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  end
  xlabel("Iterations",fontsize=label_fontsize)
  ylabel("\$\\gamma\$",fontsize=label_fontsize)
  legend(loc="best")
  # title("gamma")
  tight_layout()
  savefig( string(folder,"/gamma_all_",No_MNO,".pdf") )

  variable_convergence()
end

function variable_convergence()
  label_fontsize1 = label_fontsize
  marker_size=6
  l_width=1.2
  # stride = convert(Int, max_iters/10) +1
  stride = 6

  colors=["m","b","k","g","r","coral"]
  algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
  markers = ["x","o",">","^"]

  clf()
  figure(5,figsize=fig_size)
  subplot(3,1,1)
  plot(opt_alpha[1]*ones(max_iters),color=colors[5],linestyle=":")
  for i=2:No_MNO
    plot(opt_alpha[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  plot(alphas1[1,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width)
  for i=2:No_MNO
    plot(alphas1[i,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width)
  end
  # plot(alphas2[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[2])
  # for i=2:No_MNO
  #   plot(alphas2[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end
  plot(alphas3[1,1:max_iters],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width)
  for i=2:No_MNO
    plot(alphas3[i,1:max_iters],marker=markers[i], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width)
  end
  # plot(alphas4[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[4])
  # for i=2:No_MNO
  #   plot(alphas4[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end


  # title("alpha")
  # xlabel("Iterations",fontsize=label_fontsize1)
  ylabel("\$\\alpha\$",fontsize=label_fontsize1)
  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  # legend(bbox_to_anchor=(0.26, 0.3, 0.7, .1), loc=3,
  #        ncol=3, mode="expand", borderaxespad=0.)

  subplot(3,1,2)
  plot(opt_beta[1]*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  for i=2:No_MNO
    plot(opt_beta[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  plot(betas1[1,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width,label=algs[1])
  for i=2:No_MNO
    plot(betas1[i,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width)
  end
  # plot(betas2[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[2])
  # for i=2:No_MNO
  #   plot(betas2[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end
  plot(betas3[1,1:max_iters],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width,label=algs[3])
  for i=2:No_MNO
    plot(betas3[i,1:max_iters],marker=markers[i], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width)
  end
  # plot(betas4[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[4])
  # for i=2:No_MNO
  #   plot(betas4[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end

  # legend(loc="best")
  # title("beta")
  # xlabel("Iterations",fontsize=label_fontsize1)
  ylim(-0.05,1)
  ylabel("\$\\beta\$",fontsize=label_fontsize)
  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  legend(bbox_to_anchor=(0.13, 0.72, 0.86, .1), loc=3,
         ncol=3, mode="expand", borderaxespad=0.,  fontsize=legend_fontsize-1)

  subplot(3,1,3)
  plot(opt_gamma[1]*ones(max_iters),color=colors[5],linestyle=":")
  for i=2:No_MNO
    plot(opt_gamma[i]*ones(max_iters),color=colors[5],linestyle=":")
  end
  plot(gammas1[1,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width)
  for i=2:No_MNO
    plot(gammas1[i,1:max_iters],color=colors[2],linestyle="--",linewidth=l_width)
  end
  # plot(gammas2[1,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[2])
  # for i=2:No_MNO
  #   plot(gammas2[i,:],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end
  plot(gammas3[1,1:max_iters],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width,label=algs[3],label="Tenant 1")
  for i=2:No_MNO
    plot(gammas3[i,1:max_iters],marker=markers[i], markersize=marker_size,markevery=stride,linestyle="-",linewidth=l_width,label="Tenant $i")
  end
  # plot(gammas4[1,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width,label=algs[4])
  # for i=2:No_MNO
  #   plot(gammas4[i,:],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle=":",linewidth=l_width)
  # end

  legend(loc="best",fontsize=legend_fontsize-2)
  xlabel("Iterations",fontsize=label_fontsize1+1)
  ylabel("\$\\gamma\$",fontsize=label_fontsize1+1)
  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  savefig( string(folder,"/b_a_g_all_",No_MNO,".pdf") )

end

function read_results()
  h5open(string(folder,"results1.h5"), "r") do file
     global alphas1 = read(file, "alphas")
     global betas1 = read(file, "betas")
     global gammas1 = read(file, "gammas")
     global total_cost1 = read(file, "total_cost")

  end
  h5open(string(folder,"results2.h5"), "r") do file
     global alphas2 = read(file, "alphas")
     global betas2  = read(file, "betas")
     global gammas2 = read(file, "gammas")
     global total_cost2 = read(file, "total_cost")
  end
  h5open(string(folder,"results3.h5"), "r") do file
     global alphas3 = read(file, "alphas")
     global betas3  = read(file, "betas")
     global gammas3 = read(file, "gammas")
     global total_cost3 = read(file, "total_cost")
  end
  h5open(string(folder,"results4.h5"), "r") do file
     global alphas4 = read(file, "alphas")
     global betas4  = read(file, "betas")
     global gammas4 = read(file, "gammas")
     global total_cost4 = read(file, "total_cost")
  end

end

function plot_cost_2_scenarios()
  h5open(string(folder,"results1.h5"), "r") do file
     global total_cost1 = read(file, "total_cost")
  end

  h5open(string(folder,"results3.h5"), "r") do file
     global total_cost3 = read(file, "total_cost")
  end

  h5open(string("figs2/","results1.h5"), "r") do file
     global total_cost2_1 = read(file, "total_cost")
  end

  h5open(string("figs2/","results3.h5"), "r") do file
     global total_cost2_3 = read(file, "total_cost")
  end

  h5open(string("figs2/","results1.h5"), "r") do file
    global opt_cost1  = read(file, "opt_cost")
  end
  h5open(string("figs2/","results1.h5"), "r") do file
    global opt_solver  = read(file, "opt_solver")
  end

  marker_size=6
  l_width=1.2
  # stride = convert(Int, max_iters/10) +1
  stride = 6
  label_fontsize1 = label_fontsize

  colors=["m","b","k","g","r","coral"]
  algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
  markers = ["s","o","d"]

  figure(6,figsize=fig_size)
  subplot(2,1,1)
  plot(opt_cost1*ones(max_iters),color=colors[6],marker=markers[3], markersize=marker_size,markevery=stride,linestyle="--",label=algs[6])
  plot(opt_solver*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(total_cost2_1[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",label=algs[1])
  plot(total_cost2_3[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="-",label=algs[3])
  ylim(1,12)
  # xlabel("Iterations",fontsize=label_fontsize1)
  ylabel("Total Cost",fontsize=label_fontsize1-1)
  legend(loc="best",fontsize=legend_fontsize-1)
  title("Two Tenants Scenario")
  tight_layout(pad=1., w_pad=0.5, h_pad=0.5)

  subplot(2,1,2)
  plot(opt_cost*ones(max_iters),color=colors[5],linestyle=":",label=algs[5])
  plot(total_cost1[1:max_iters],color=colors[2],marker=markers[1], markersize=marker_size,markevery=stride,linestyle="-",label=algs[1])
  plot(total_cost3[1:max_iters],color=colors[4],marker=markers[2], markersize=marker_size,markevery=stride,linestyle="-",label=algs[3])

  xlabel("Iterations",fontsize=label_fontsize1)
  ylabel("Total Cost",fontsize=label_fontsize1-1)
  legend(loc=1,fontsize=legend_fontsize-1)
  ylim(1,12)
  title("Four Tenants Scenario")
  tight_layout(pad=1., w_pad=0.5, h_pad=0.5)

  savefig( string(folder,"wrap_cost_Convergence_all_",No_MNO,".pdf") )
end
