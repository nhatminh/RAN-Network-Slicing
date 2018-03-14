# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:19:14 2017

@author: Minh
"""

import h5py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

N_Sims= 300

#### DECENTRALIZE PLOT ####
def plot_decentralize():
    filename = 'figs4\MC{0}.h5'.format(N_Sims)
    f = h5py.File(filename, 'r')
    # Get the data
    for name in f:
        print(name)
        
    Cost_MC = f['Cost_MC'][:,:,:]
    Converged_Iters_MC  = f['Converged_Iters_MC'][:,:]
    print(Cost_MC.shape)
    print(Converged_Iters_MC.shape)
    
    cost_cols  =['Utilization','Backhaul','Total']
    iters_cols =['PBCD','JP-ADMM']
    df_Cost_alg1 = pd.DataFrame(Cost_MC[:,:,0], columns=cost_cols)
    df_Cost_alg3 = pd.DataFrame(Cost_MC[:,:,2], columns=cost_cols)
    df_iters = pd.DataFrame(Converged_Iters_MC[:,[0,2]], columns=iters_cols)
    
    print(df_Cost_alg1.head())
    #print(df_Cost_alg3.head())
    print(df_iters.head())
    
    df_Cost_alg1 = pd.melt(df_Cost_alg1,var_name='Cost')
    df_Cost_alg3 = pd.melt(df_Cost_alg3,var_name='Cost')
    df_iters = pd.melt(df_iters,var_name='Alg',value_name='Iterations')
    print(df_Cost_alg1.head())
    print(df_iters.head())
    
    df_Cost_alg1['Alg'] = 'PBCD'
    df_Cost_alg3['Alg'] = 'JP-ADMM'
    merge_df = pd.concat([df_Cost_alg1,df_Cost_alg3])
    print(merge_df.head())
    
    plt.subplot(1,3,1)
    sns.violinplot(x="Cost",y="value", data=df_Cost_alg1, palette="muted")
    plt.subplot(1,3,2)
    sns.swarmplot(x="Cost",y="value", data=df_Cost_alg1)
    plt.subplot(1,3,3)
    sns.boxplot(x="Cost",y="value", data=df_Cost_alg1)
    plt.show()
    
    plt.subplot(1,3,1)
    sns.violinplot(x="Cost",y="value",hue="Alg", data=merge_df, palette="muted")
    plt.subplot(1,3,2)
    sns.swarmplot(x="Cost",y="value",hue="Alg", data=merge_df)
    plt.subplot(1,3,3)
    sns.boxplot(x="Cost",y="value",hue="Alg", data=merge_df)
    plt.show()
    
    
    plt.subplot(1,2,1)
    sns.violinplot(x="Alg",y="Iterations", data=df_iters, palette="muted")
    plt.subplot(1,2,2)
    sns.swarmplot(x="Alg",y="Iterations", data=df_iters)
    plt.show()
    
    plt.subplot(1,2,1)
    sns.boxplot(x="Alg",y="Iterations", data=df_iters)
    plt.subplot(1,2,2)
    sns.barplot(x="Alg",y="Iterations", data=df_iters)
    plt.show()
    
  
### CENTRALIZE PLOT ####
def plot_centralize1():
    filename = 'figs4\MC{0}-Cen.h5'.format(N_Sims)
    f = h5py.File(filename, 'r')
    # Get the data
    for name in f:
        print(name)
        
    Cost_MC_Cen = f['Cost_MC_Cen'][:,:]
    print(Cost_MC_Cen.shape)
    
    cost_cols  =['U_Cost','B_Cost','T_Cost']
    df_Cost = pd.DataFrame(Cost_MC_Cen, columns=cost_cols)
    print(df_Cost.head())
    
    df_Cost = pd.melt(df_Cost,var_name='Cost')
    print(df_Cost.head())
    
    plt.subplot(1,3,1)
    sns.violinplot(x="Cost",y="value", data=df_Cost, palette="muted")
    plt.subplot(1,3,2)
    sns.swarmplot(x="Cost",y="value", data=df_Cost)
    plt.subplot(1,3,3)
    sns.boxplot(x="Cost",y="value", data=df_Cost)
    plt.show()
    

### CENTRALIZE PLOT ####
def plot_centralize2():
    w_tradeoff=[0.01,0.1,1]
    
    filename = 'figs4\MC{0}-Cen.h5'.format(N_Sims)
    f = h5py.File(filename, 'r')
    # Get the data
    for name in f:
        print(name)
        
    for w in range(len(w_tradeoff)):  
        Cost_MC_Cen = f['Cost_MC_Cen{0}'.format(w+1)][:,:]
        print(Cost_MC_Cen.shape)
        
        cost_cols  =['U_Cost','B_Cost','T_Cost']
        df_Cost = pd.DataFrame(Cost_MC_Cen, columns=cost_cols)
        print(df_Cost.head())
        
        df_Cost = pd.melt(df_Cost,var_name='Cost')
        print(df_Cost.head())
        
        plt.subplot(1,3,1)
        sns.violinplot(x="Cost",y="value", data=df_Cost, palette="muted")
        plt.subplot(1,3,2)
        sns.swarmplot(x="Cost",y="value", data=df_Cost)
        plt.subplot(1,3,3)
        sns.boxplot(x="Cost",y="value", data=df_Cost)
        plt.show()
        
### CENTRALIZE PLOT ####
def plot_final():
    cost_cols  =['Aggregated BS Utilization','Backhaul Link Cost','Total Cost']
        
    df_u, df_b, df_t, df_iters = read_files()
#    print(df_u.head())
    df_u = df_u.loc[df_u[cost_cols[0]] <= 4]
    df_b = df_b.loc[df_b[cost_cols[1]] <= 30]
    df_t = df_t.loc[df_t[cost_cols[2]] <= 15]
    print(df_u.head())
    
    plt.figure(1,figsize=(7.,5.1))
    sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.7)
    
    plt.subplot(1,2,1)
    sns.boxplot(x='$\omega$',y=cost_cols[0], data=df_u,width=0.55 )
    plt.subplot(1,2,2)
    sns.boxplot(x='$\omega$',y=cost_cols[1], data=df_b,width=0.55 )
#    plt.subplot(1,3,3)
#    sns.boxplot(x='$\omega$',y=cost_cols[2], data=df_t,width=0.55 )
#    sns.despine()
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("figs4\Cost_Stat_300.pdf")
#    plt.savefig("figs4\Cost_Stat_100.pdf")
    plt.show()
    
    plt.figure(2,figsize=(7.,5.1))
    sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.7)
    sns.swarmplot(x="Algorithm",y="Iterations", data=df_iters)
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.xlabel('')
    plt.savefig("figs4\Iters_Stat_300.pdf")
#    plt.savefig("figs4\Iters_Stat_100.pdf")
    plt.show()


def read_files(): 
#    filenames = [ 'figs4\MC100-Cen.h5', 'figs4\MC200-Cen.h5', 'figs4\MC200-Cen1.h5']
    filenames = [ 'figs4\MC100-Cen.h5', 'figs4\MC200-Cen.h5']
#    filenames = [ 'figs4\MC{0}-Cen.h5'.format(N_Sims)]
    w_tradeoff=['0.01','0.1','1']
    
    df_iters    = []
    df_util     = []
    df_backhaul = []
    df_total    = []
    
    for file in filenames:
        f = h5py.File(file, 'r')
        # Get the data
        for name in f:
            print(name)
            
        Converged_Iters_MC  = f['Converged_Iters_MC'][:,:]
        print(Converged_Iters_MC.shape)
        
        cost_cols  =['Aggregated BS Utilization','Backhaul Link Cost','Total Cost']
        iters_cols =['PBCD','JP-ADMM']
        
        df_iter = pd.DataFrame(Converged_Iters_MC[:,[0,2]], columns=iters_cols)    
        df_iter = pd.melt(df_iter,var_name='Algorithm',value_name='Iterations')
        df_iters.append(df_iter)
     
        for w in range(len(w_tradeoff)):
            Cost_MC_Cen = f['Cost_MC_Cen{}'.format(w+1)][:,:]
            df_Cost = pd.DataFrame(Cost_MC_Cen, columns=cost_cols)
            df_Cost['$\omega$'] = w_tradeoff[w]
            print(df_Cost.tail())
#            df_Cost = pd.melt(df_Cost,var_name='Cost',value_name='Value')

            df_util.append(df_Cost[[cost_cols[0],'$\omega$']])
            df_backhaul.append(df_Cost[[cost_cols[1],'$\omega$']])
            df_total.append(df_Cost[[cost_cols[2],'$\omega$']])
        

#
#        Cost_MC_Cen2 = f['Cost_MC_Cen2'][:,:]
#        df_Cost2 = pd.DataFrame(Cost_MC_Cen2, columns=cost_cols)
#        df_Cost2 = pd.melt(df_Cost2,var_name='Cost',value_name='Value')
#        df_w2.append(df_Cost2)
#        
#        Cost_MC_Cen3 = f['Cost_MC_Cen3'][:,:]
#        df_Cost3 = pd.DataFrame(Cost_MC_Cen3, columns=cost_cols)
#        df_Cost3 = pd.melt(df_Cost3,var_name='Cost',value_name='Value')
#        df_w3.append(df_Cost3)

    return pd.concat(df_util), pd.concat(df_backhaul), pd.concat(df_total), pd.concat(df_iters) 
        
        
#plot_centralize2()
plot_final()