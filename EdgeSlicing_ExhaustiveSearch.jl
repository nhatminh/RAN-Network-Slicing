using Convex
using ECOS
using Gurobi
using SCS
using Ipopt
using PyPlot
using HDF5

global REUSED_INPUT = true
global READ_A_G = false
global READ_A_B= true
F= 5*10^4
S= 220
B= 500
L= 100*1024*8  #File size
No_MNO=2
max_iters = 55
ALG_MODE = 4 #1, 2, 3, 4
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
   end
else
  include("Setting2.jl")
  h5open(input_filename, "w") do file
      write(file, "lamb_n", lamb_n)
      write(file, "w_ulti", w_ulti)
      write(file, "w_lamb_hit_n", w_lamb_hit_n)
  end
end

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

rho = 0.05
rho1 = 10.
# lamb_n = 0
# w_ulti = [0.9, 0.9 ]
# w_lamb_hit = [1., 1.]

opt_alpha = zeros(n)
opt_beta  = zeros(n)
opt_gamma = zeros(n)


function exhaustive_search(interval)
  s_alphas = zeros(No_MNO)
  s_betas  = zeros(No_MNO)
  s_gammas = zeros(No_MNO)
  global opt_cost=1000

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
              if (total_cost < opt_cost)
                opt_alpha[:] = s_alphas[:]
                opt_beta[:]  = s_betas[:]
                opt_gamma[:] = s_gammas[:]
                opt_cost = total_cost
              end
            end
          end
      end
    end
  end
  println("alpha: ",opt_alpha)
  println("beta: " ,opt_beta)
  println("gamma: ",opt_gamma)
  println("Opt_cost: ",opt_cost)
  alpha_beta_gamma_relationship()
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
  if (sum(alphas) !=1)  || (sum(betas) !=1) || (sum(gammas) !=1)
    fes_Flag= false
  end

  for i = 1:n
    if (alphas[i] == 0) || (alphas[i] < w_ulti[i]*betas[i]/(1-epsilon))
      fes_Flag= false
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
function alpha_beta_gamma_relationship()
  lamb_miss = max(0,lamb_n - w_lamb_hit_n .*opt_beta)
  println( "miss cost: ", lamb_miss./(opt_gamma[:]*B - lamb_miss))
  println( "utilization: ", w_ulti .*(opt_beta./opt_alpha))
end


using Ipopt
using PyPlot
using JuMP

interval=0.02
numb_slices = convert(Int,1/interval)
a_g_cost = zeros(numb_slices,numb_slices)
a_b_cost = zeros(numb_slices,numb_slices)
fig_size = (7.,5.1)
label_fontsize = 16
legend_fontsize = label_fontsize - 4

function centralized_solver1(alpha, beta)
  prob = Model(solver=IpoptSolver(tol=1e-9, max_iter=10000, print_level =0))
  @variable(prob, gamma[1:No_MNO] >= 0)
  @variable(prob, lamb_miss[1:No_MNO]>=0)

  @NLobjective(prob, Min,  sum(w_ulti[i]*beta[i]/alpha[i] for i=1:No_MNO) + w_tradeoff*( sum(lamb_miss[i]/(gamma[i]*B-lamb_miss[i]) for i=1:No_MNO) ) )

  @constraint(prob, sum(gamma[i] for i =1:No_MNO) == 1 )

  for i = 1:No_MNO
    @constraint(prob,  w_ulti[i]*beta[i]  - (1 - epsilon)*alpha[i] <= 0 )
    @constraint(prob, lamb_miss[i] +1e-8 <= gamma[i]*B )
    @constraint(prob, lamb_miss[i] == lamb_n[i] - w_lamb_hit_n[i] * beta[i])
  end

  status = solve(prob)
  # println("Solve Status: ",status)
  g_vals = getvalue(gamma)
  if(status==:Optimal)
    return  g_vals, min(10,cal_cost_MVNOs(alpha,beta,g_vals))
  else
    return g_vals, 10
  end
end

function save_a_b()
  h5open("figs2/a_b_map.h5", "w") do file
      write(file, "cost", a_b_cost)
      write(file, "numb_slices", numb_slices)
  end
end

function read_a_b()
  h5open("figs2/a_b_map.h5", "r") do file
      global a_b_cost = read(file, "cost")
      global numb_slices = read(file, "numb_slices")
  end
end

function read_a_b_decentralized()
  h5open("figs2/results1.h5", "r") do file
     global alphas2 = read(file, "alphas")
     global betas2 = read(file, "betas")
     global total_cost2 = read(file, "total_cost")
  end
  h5open("figs2/results3.h5", "r") do file
     global alphas4 = read(file, "alphas")
     global betas4 = read(file, "betas")
     global total_cost4 = read(file, "total_cost")
  end
  alg_interval = 5
  global converge_k2 = 1
  for i=1:size(alphas2[1,:])[1]
    if((abs(alphas2[1,i] - alphas2[1,i+1]) + abs(alphas2[1,i]- alphas2[1,i+1])) <1e-4 )
      break
    end
    converge_k2 += 1
  end
  println(converge_k2)
  converge_k2 = max_iters

  global converge_k4 = 1
  for i=1:size(alphas4[1,:])[1]
    if((abs(alphas4[1,i] - alphas4[1,i+1]) + abs(alphas4[1,i]- alphas4[1,i+1])) <1e-4 )
      break
    end
    converge_k4 += 1
  end
  println(converge_k4)
  converge_k4 = max_iters
end

function plot_a_b_map()
  read_a_b_decentralized()

  if(READ_A_B)
    read_a_b()
  else
    for i = 1:numb_slices
      alpha1= interval*i
      alpha2= 1-alpha1
      for j = 1:numb_slices
        beta1= interval*j
        beta2= 1- beta1
        _,a_b_cost[i,j] = centralized_solver1([alpha1,alpha2],[beta1,beta2])
      end
    end
    save_a_b()
  end

  for i = 1:numb_slices
    for j = 1:numb_slices
      a_b_cost[i,j] = min(a_b_cost[i,j],4)
    end
  end

  x_start2 = betas2[1,1]
  y_start2 = alphas2[1,1]
  x_end2 = betas2[1,converge_k2]
  y_end2 = alphas2[1,converge_k2]
  x_start4 = betas4[1,1]
  y_start4 = alphas4[1,1]
  x_end4 = betas4[1,converge_k4]
  y_end4 = alphas4[1,converge_k4]

  alpha = linspace(0, 1, numb_slices)
  beta = linspace(0, 1, numb_slices)

  xgrid = repmat(alpha',numb_slices,1)
  ygrid = repmat(beta,1,numb_slices)

  figure("pyplot_surfaceplot1",figsize=fig_size)
  # ax = fig[:add_subplot](2,1,1, projection = "3d")
  #cmap=ColorMap("YlOrRd")
  # ax[:plot_surface](xgrid, ygrid, a_g_cost, rstride=3,edgecolors="k", cstride=3, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
  # ax[:contourf](xgrid, ygrid, a_g_cost, zdir='z', offset=-100, cmap=ColorMap("coolwarm"))
  # scatter3D(gammas[1,1:converge_k],alphas[1,1:converge_k],a_g_cost[1:converge_k],s=35,color="b",alpha=0.4)

  plot_surface(xgrid, ygrid, a_b_cost, rstride=3,edgecolors="k", cstride=3, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
  xlabel("\$\\beta\$",fontsize=label_fontsize)
  ylabel("\$\\alpha\$",fontsize=label_fontsize)
  tight_layout()
  savefig("figs2/a_b_map_3D.pdf")


  fig = figure("pyplot_surfaceplot3",fig_size)
  # cp = contour(xgrid, ygrid, a_b_cost,cmap="terrain")
  cp = contour(xgrid, ygrid, a_b_cost,cmap="coolwarm")
  clabel(cp, inline=1, fontsize=11)
  # scatter(gammas4[1,1:converge_k4],alphas4[1,1:converge_k4],s=45,marker="o",color="r",alpha=0.4,label="Decentralized")
  # scatter(gammas2[1,1:converge_k2],alphas2[1,1:converge_k2],s=45,marker="^",color="b",alpha=0.4,label="Centralized")
  plot(betas4[1,1:2:converge_k4],alphas4[1,1:2:converge_k4],linestyle="-",marker="o",markersize=8, color="darkgreen",alpha=0.6,label="JP-ADMM")
  plot(betas2[1,1:2:converge_k2],alphas2[1,1:2:converge_k2],linestyle="--",marker="^",markersize=8,color="blueviolet",alpha=0.5,label="PBCD")

  annotate("Starting Point",
            xy=[x_start2;y_start2],
            xytext=[x_start2+0.07;y_start2+0.07],
            xycoords="data",
            size=14,
            arrowprops=Dict("arrowstyle"=>"fancy",
                "facecolor"=>"orangered",
                "connectionstyle"=>"angle3,angleA=0,angleB=-90"))

  annotate("Optimal",
            xy=[x_end2;y_end2],
            xytext=[x_end2-0.19;y_end2-0.08],
            xycoords="data",
            size=14,
            arrowprops=Dict("arrowstyle"=>"fancy",
                "facecolor"=>"red",
                "edgecolor"=>"black",
                "connectionstyle"=>"angle3,angleA=0,angleB=+45"))


  legend(loc=2,fontsize=legend_fontsize+2)
  xlabel("\$\\beta\$",fontsize=label_fontsize+2)
  ylabel("\$\\alpha\$",fontsize=label_fontsize+2)
  tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
  savefig("figs2/a_b_map_contour2.pdf")
end

plot_a_b_map()
