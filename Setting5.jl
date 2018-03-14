using Distributions

InP=[250,250];

# User generation
No_MNO=5
N_UE_s=[20,19,19,18,20]
ReRate_MNO=[15,16,16,16,15]
# No_MNO=5
# N_UE_s=[20,18,20,15,25]
# ReRate_MNO=[15,12,18,10,16]
# F=5*10^4
# S=220
# B=1000
# L= 100*1024*8  #File size

N_UE=sum(N_UE_s,1)[1]
UE = zeros(2,N_UE)
uniform_loc = Uniform(1,500)
UE[1,1:N_UE]=rand(uniform_loc,N_UE)
UE[2,1:N_UE]=rand(uniform_loc,N_UE)

# BS and UE distance
d=zeros(N_UE,1)
for i = 1:N_UE
  d[i,1]=norm(InP-UE[:,i])
end
d;

for i = 1:N_UE
  if(d[i,1]==0)
      d[i,1]=d[i,1]+randi([1,500]);
  end
end

# pathloss calculation
PL = zeros(N_UE);
for i = 1:N_UE
  PL[i]=34 + 40*log10(d[i,1]/1000);
end

# channel gain calculation
normal_distribution = Normal(0,8)
H = zeros(N_UE)
for i = 1:N_UE
    H[i]= PL[i]+ rand(normal_distribution)
end

# Pr
Pr = zeros(N_UE)
for i = 1:N_UE
    Pr[i]= 49 - H[i]
end

# C_u
gxP = 10.^(Pr/10)
SNR=zeros(N_UE)
C_u=zeros(N_UE)

for i = 1:N_UE
   SNR[i]= gxP[i]/(20*10^(6-10.4))
   C_u[i]=20*10^6*log2(1+SNR[i])
end

# matrix lambda_u,f
Zipf_parameter =rand(N_UE)
prob_u_f_temp=zeros(N_UE, F)

for i =1:N_UE
    for j=1:F
        prob_u_f_temp[i,j]=1/(j^Zipf_parameter[i])
    end
end

prob_sum=1./(sum(prob_u_f_temp,2))
reprob_sum=repmat(prob_sum, 1,F)
prob_u_f= prob_u_f_temp .* reprob_sum
rate_vt=[ReRate_MNO[1]*ones(N_UE_s[1],1); ReRate_MNO[2]*ones(N_UE_s[2],1); ReRate_MNO[3]*ones(N_UE_s[3],1);
          ReRate_MNO[4]*ones(N_UE_s[4],1); ReRate_MNO[5]*ones(N_UE_s[5],1)]
# rate_vt=[15*ones(N_UE_s[1],1); 12*ones(N_UE_s[2],1); 18*ones(N_UE_s[3],1); 10*ones(N_UE_s[4],1); 16*ones(N_UE_s[5],1)]
rate_matrix=repmat(rate_vt, 1,F)
lamb_u_f=rate_matrix.*prob_u_f
lamb_u=sum(lamb_u_f,2)
lamb_n = ReRate_MNO.*N_UE_s

# index bat dua user cua tung mvno
index=[1,21,40,59,77,97]
constant_hf=zeros(No_MNO,F)
sum_n=zeros(No_MNO,F)

for n=1:No_MNO
  sum_n[n,:]=sum(lamb_u_f[index[n]:index[n]+N_UE_s[n]-1,:],1)
  constant_hf[n,:]=sum_n[n,:] /lamb_n[n]*S
end

w_lamb_hit_u=zeros(N_UE,1)

for i=1:N_UE
    if i<21
        m=1
    elseif i<40
        m=2
    elseif i<59
        m=3
    elseif i<77
        m=4
    elseif i<97
        m=5
    elseif i<115
        m=6
    else
        m=7
    end
    w_lamb_hit_u[i]= sum(sum_n[m,:].*constant_hf[m,:],2)[1]
end

w_mu_u = C_u/L
w_ulti = zeros(No_MNO)
w_lamb_hit_n = zeros(No_MNO)
w_ulti_u = w_lamb_hit_u./w_mu_u
for n=1:No_MNO
  w_ulti[n] = sum(w_ulti_u[index[n]:index[n]+N_UE_s[n]-1])
  w_lamb_hit_n[n] = sum(w_lamb_hit_u[index[n]:index[n]+N_UE_s[n]-1])
end

# println(w_ulti)
# # println(w_lamb_hit_u)
# println(w_lamb_hit_n)
