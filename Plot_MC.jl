fig_size = (7.,5.1)
label_fontsize = 17
legend_fontsize = label_fontsize - 4
patterns = ["","."]

folder = string("figs$No_MNO/")
using PyPlot
# using Seaborn

# function plt_MC_Hist()
#   # figure(1)
#   # plt[:hist](Cost_MC[1,1,:], bins=20, normed=true, alpha=0.5)
#   # plt[:hist](Cost_MC[1,2,:], bins=20, normed=true, alpha=0.5)
#   # plt[:hist](Cost_MC[1,3,:], bins=20, normed=true, alpha=0.5)
#
# #   figure(2)
# #   hist(Cost_MC[3,1,:], bins=20, normed=true, color='r', alpha=0.3)
# #   hist(Cost_MC[3,2,:], bins=20, normed=true, color='b', alpha=0.3)
# #   hist(Cost_MC[3,3,:], bins=20, normed=true, color='g', alpha=0.3)
# #
#   figure(3)
#   plt[:hist](Converged_Iters_MC[1,:], bins=20, normed=true, alpha=0.5)
#   figure(4)
#   plt[:hist](Converged_Iters_MC[3,:], bins=20, normed=true, alpha=0.5)
# end
#
# function plt_Stats()
#   figure(5)
#   Seaborn.boxplot(Converged_Iters_MC[1,:])
#   figure(6)
#   Seaborn.boxplot(Converged_Iters_MC[3,:])
#   figure(7)
#   Seaborn.boxplot(Converged_Iters_MC)
# end

function read_MC()
  # h5open(string(folder,"MC$N_SIMs-Cen-eps1.h5"), "r") do file
  h5open(string(folder,"MC$N_SIMs-Cen-eps.h5"), "r") do file
     global Cost_MC_Cen3_e1 = read(file, "Cost_MC_Cen3_e1")
     global Cost_MC_Cen3_e2 = read(file, "Cost_MC_Cen3_e2")
  end
end

function plt_epsilon()
  fig_size = (7.,5.1)
  label_fontsize = 15
  legend_fontsize = label_fontsize - 3

  figure(1,figsize=fig_size)
  set_style("whitegrid")
  subplot(1,2,1)
  cost1 = zeros(No_MNO)
  cost2 = zeros(No_MNO)
  cost3 = zeros(No_MNO)
  idx= 1:4
  # label_operators = ["Operator 1", "Operator 2", "Operator 3", "Operator 4"]
  label_operators = ["1", "2", "3", "4"]
  for i = 1:No_MNO
    # cost1[i] = median(Cost_MC_Cen3_e1[1,i,:])
    # cost2[i] = median(Cost_MC_Cen3_e1[2,i,:])
    cost1[i] = mean(Cost_MC_Cen3_e1[1,i,:])
    cost2[i] = mean(Cost_MC_Cen3_e1[2,i,:])
    cost3[i] = cost1[i] + cost2[i]
    # cost3[i] = cost1[i] + 0.1*cost2[i]
  end
  b1 = bar(idx,cost1)
  b2 = bar(idx,cost2, bottom=cost1)
  xticks(idx,label_operators)
  legend((b1[1],b2[1]), ("Utilization","Backhaul"),loc=1,fontsize=legend_fontsize)
  ylabel("Cost",fontsize=label_fontsize)
  xlabel("Tenant",fontsize=label_fontsize)
  title("\$\\epsilon = 0.001\$")
  # ylim(0,1.2)
  # ylim(0,1.6)
  ylim(0,2.2)
  println("Total_cost:",cost3 )
  println("Total_cost2:",sum(cost3) )

  subplot(1,2,2)
  for i = 1:No_MNO
    # cost1[i] = median(Cost_MC_Cen3_e2[1,i,:])
    # cost2[i] = median(Cost_MC_Cen3_e2[2,i,:])
    cost1[i] = mean(Cost_MC_Cen3_e2[1,i,:])
    cost2[i] = mean(Cost_MC_Cen3_e2[2,i,:])
    cost3[i] = cost1[i] + cost2[i]
    # cost3[i] = cost1[i] + 0.1*cost2[i]
  end
  b1 = bar(idx,cost1)
  b2 = bar(idx,cost2, bottom=cost1)
  xticks(idx,label_operators)
  legend((b1[1],b2[1]), ("Utilization","Backhaul"),loc=1,fontsize=legend_fontsize)
  ylabel("Cost",fontsize=label_fontsize)
  xlabel("Tenant",fontsize=label_fontsize)
  title("\$\\epsilon = 0.1\$")
  # ylim(0,1.2)
  # ylim(0,1.6)
  ylim(0,2.2)
  println("Total_cost:",cost3 )
  println("Total_cost2:",sum(cost3) )

  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  # savefig( string(folder,"cost_comparison_epsilon_",No_MNO,".pdf") )
  savefig( string(folder,"cost_comparison_epsilon_",No_MNO,"_1.pdf") )

end
