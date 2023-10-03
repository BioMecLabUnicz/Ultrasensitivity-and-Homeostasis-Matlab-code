% The script generates the results shown in Figures 4 and 5. 
% In the script there is a flag (tag_case) that allows generating
% the results of Fig 4 and Fig 5:
% tag_case=1 -> Fig. 4
% tag_case=2 -> Fig. 5
% For each Figure, the flag tag_nh allows to get the results for different
% values of the exponents n_{Hog1} and n_{Fps1}: 
% tag_nh=1; n_{Hog1}=n_{Fps1}=2, plot the panels of first column of Figs. 4C-D or  Figs. 5C-D
% tag_nh=2; n_{Hog1}=n_{Fps1}=4, plot the panels of second column of Figs. 4C-D or  Figs. 5C-D
% tag_nh=3; n_{Hog1}=n_{Fps1}=8, plot the panels of third column of panels Figs. 4C-D or  Figs. 5C-D
% For generating the steady state curves you need to simulate the ode model
% (described by the function yeast_model.m)
% without perturbation (the perturbation time global_t0=0 and the relative
% amplitude global_PIe_d=0), and varying the input signal ka (by the vector
% vct_Ks); then by computing the steady state variables for each input, you
% can obtain the input-output curves.


clear all
close all
clc

tag_case=1; % (K_hog1 and K_fps1=0.034 microM)
% tag_case=2; % (K_hog1 and K_fps1=0.017 microM)

tag_nh=1; % (n_h=n_{Hog1}=n_{Fps1}=2)
% tag_nh=2; % (n_h=n_{Hog1}=n_{Fps1}=4)
% tag_nh=3; % (n_h=n_{Hog1}=n_{Fps1}=8)

%define global paramters
global ka ks krs km krm khs khm H_0 Rs_0 Rm_0 RCtot_0 HPTtot_0
global krc kc krt kt khc 


%%% sln1 system parameters
c_dim=1e3; % constant used for converting from microM to milliM
c_time=60; % constant used for converting from sec^{-1} to min^{-1}
kc=160*c_dim*c_time;
krc=0*c_dim*c_time;

kt=20.7*c_dim*c_time;
krt=29.5*c_dim*c_time;

ks=66.67*c_dim*c_time;
km=1*c_dim*c_time;
 
%reverse phosphorylation rates 
krs=0;
krm=0.08*c_dim*c_time;

% hydrolysis rate/autodephosphorylation rate

khc=0.05*c_time;
khs=0.05*c_time;
khm=0.05*c_time;

%total protein concentrations in all layers
H_0=0.25/c_dim;
RCtot_0=0.25/c_dim;
HPTtot_0=1.5/c_dim;
Rm_0=1.5/c_dim;
Rs_0=1.5/c_dim;


% define global parameters for mapk hog1 system
global S_0 a1 a2 E2tot_0 K1 K2 k1 k2 d1 d2 Hog1_0;

% define global parameters for glycerol regulation
global  n_hog_fps1   n_hog_glyc K_hog_glyc  K_hog_fps1 

S_0=1.5*1e-3;         
a1=100*1e6*c_time;  
a2=100*1e6*c_time; 
K1=0.010*1e-3;    
K2=0.010*1e-3;    
k1=1*c_time;
k2=1*c_time;

d1=K1*a1-k1;
d2=K2*a2-k2;


% perturbation  at time 0, with amplitude PIe_d
global global_t0 global_PIe_d

% time_0
global_t0=0; % can be a vector where each element represents the time of perturbation

% amplitude of perturbation
PIe_d=1.86*0.4;
global_PIe_d= PIe_d; % can be a vector where each element represents the amplitude of perturbation

% global control and model parameters
global global_V_b global_V_e global_V_PIt0 global_V_0;
global global_PIi_0 global_PIe_0 global_glyc_0;
global global_kp1
global global_k_hog global_kfps1

%nominal values for the model parameters - see Table S1
global_kp1=100; % water perm. constant              
global_V_b=0.368; % non osmotic volume[0.31 0.46]
global_V_e=4.79e3; % external volume [0.5e3 5e3]
global_PIi_0=0.65; % initial internal pressure [0.6 0.7]
global_PIe_0=0.240; % initial external pressure [0.24 0.25]
global_glyc_0=2e-4; % initial glycerol [1.1e-4 5e-4]
PIt_0=global_PIi_0-global_PIe_0; % initial turgor pressure
global_V_PIt0=0.85; % volume when the turgor pressure is zero [0.75 0.95]
global_V_0=1;

global_kfps1=.2;
E2tot_0=.425/c_dim;  

if  tag_case==1
    K_hog_fps1=.08*E2tot_0;
    K_hog_glyc=.08*E2tot_0;
    global_k_hog=0.021;
elseif  tag_case==2
    K_hog_fps1=.04*E2tot_0;  
    K_hog_glyc=.04*E2tot_0;
    global_k_hog=0.019;
else
    K_hog_fps1=.04*E2tot_0;  
    K_hog_glyc=.04*E2tot_0;
    global_k_hog=0.019;
end


l1_w=165; h1_w=165;
if tag_nh==1
    n_hog_fps1=2; 
    n_hog_glyc=2;
    
    if tag_case==1
       vct_Ks=[0.7500    0.8500    1.0000    1.1250    1.2500];
        %vct_K8_n2_fin=[0.7500    0.8500    1.0000    1.1250    1.2500];

    elseif tag_case==2        
        vct_Ks=[0.8250    0.9750    1.1000    1.2500    1.4000];
        %vct_K4_n2_fin=[0.8250    0.9750    1.1000    1.2500    1.4000];

    end
elseif  tag_nh==2
    n_hog_fps1=4; 
    n_hog_glyc=4;
    
    if tag_case==1
       vct_Ks=[0.7250    0.8250    0.9750    1.0750    1.2000];
       %vct_K8_n4_fin=[0.7250    0.8250    0.9750    1.0750    1.2000];
 
    elseif tag_case==2        
        vct_Ks=[0.7750    0.9250    1.0500    1.2250    1.3250];
        %vct_K4_n4_fin=[0.7750    0.9250    1.0500    1.2250    1.3250 ];

    end
    
elseif  tag_nh==3
    n_hog_fps1=8; 
    n_hog_glyc=8;
    
    if tag_case==1
       vct_Ks=[0.7500    0.8250    0.9250    1.0250    1.1000];
        %    vct_K8_n8_fin=[0.7500    0.8250    0.9250    1.0250    1.1000 ];


    elseif tag_case==2        
        vct_Ks=[0.8500    0.9250    1.0500    1.1750    1.2750];
        %    vct_K4_n8_fin=[0.8500    0.9250    1.0500    1.1750    1.2750];

    end


else
    n_hog_fps1=2; 
    n_hog_glyc=2;
    
    if tag_case==1
       vct_Ks=[0.7500    0.8500    1.0000    1.1250    1.2500];
        %vct_K8_n2_fin=[0.7500    0.8500    1.0000    1.1250    1.2500];

    elseif tag_case==2        
        vct_Ks=[0.8250    0.9750    1.1000    1.2500    1.4000];
        %vct_K4_n2_fin=[0.8250    0.9750    1.1000    1.2500    1.4000];

    end
end
num_fig=10*(tag_case-1)+2*(tag_nh-1)+1;





%%%%INITIAL CONDITIONS %%%%
HKp0=H_0;
RCp0=RCtot_0;
HPTp0=HPTtot_0;
R1p0=Rs_0;
R2p0=Rm_0;
A0=0.1*1e-3;
C10=0;
C20=0;
Af0=0*1e-3;
C1f0=0;
C2f0=0;
Sf_tot_0=1.5*1e-3;

V0=global_V_0;
glyc0=global_glyc_0;
glyc_ex0=global_V_e*global_glyc_0/(global_V_0-global_V_b);

Hog1_0=A0;
x0=[HKp0;RCp0;HPTp0;R1p0;R2p0;A0;C10;C20;V0;glyc0;glyc_ex0;H_0;...
    RCtot_0;HPTtot_0;Rs_0;Rm_0;S_0;E2tot_0];

% data time
tspan4data=[-20,-10,-5,-2.0000, -1.3333,-0.6667,-0:.01:40,40.1:0.1:60,61:80];
% simulation time
tspan=[-500:-450, -400:50:-50, tspan4data, 100:20:200];
options=[];

count_k=0; % count variable for the different values of vct_Ks

for idx_a=vct_Ks%

    count_k=count_k+1;
    if count_k==1 
        type_col_line='-g';
    elseif  count_k==2
        type_col_line='-m';
    elseif  count_k==3
        type_col_line='-k';

    elseif  count_k==4
        type_col_line='-c';
    elseif  count_k==5
        type_col_line='-b';

    elseif  count_k==6
        type_col_line='-.k';
       
    else
        type_col_line='--g';
    end
    ka=idx_a*c_time;
    [time, Y]=ode15s('yeast_model', tspan, x0, options);

    if count_k==1
            col_y='b';
        elseif count_k==2
            col_y='k';
        elseif count_k==3
            col_y='r';
    end
    
    t_in_v=find(time>=-5 & time<=0); % time before perturbation
       
    t_in=t_in_v(1);
    V1_ss=Y(t_in,9);
%     V1_ss_v(count_k)=V1_ss;
    Rs_tot=Y(:,15);
    R1p=Y(:,4);
    C1=Y(:,7);
    R1=Rs_tot-R1p-C1;
    A=Y(:,6);    
    C2=Y(:,8);
    S_tot=Y(:,17);    
    I=S_tot-A-C1-C2;
 
    Af=A;
    Ah=A;

    u_hog=global_k_hog*Ah.^n_hog_glyc./(Ah.^n_hog_glyc+K_hog_glyc.^n_hog_glyc);
    u_Fps1=global_kfps1*(K_hog_fps1^n_hog_fps1)./(Af.^n_hog_fps1+K_hog_fps1^n_hog_fps1);
      
    f_Ah=Ah.^n_hog_glyc./(Ah.^n_hog_glyc+K_hog_glyc.^n_hog_glyc);
    f_Af=(K_hog_fps1^n_hog_fps1)./(Af.^n_hog_fps1+K_hog_fps1^n_hog_fps1);
    
    [PIi, PIe, PIt,  Vdot, Vratio]=compute_pi_vdot_vratio(time,Y);
       
    l_w=1;    
    f1=figure (num_fig);
    hold on
    grid on
    V=Y(:,9);
    idx_time_less0=find(time<=0);
    idx_time_0=idx_time_less0(end);
    V0_max=V(idx_time_0);
    plot(time,V/V0_max,type_col_line,'LineWidth',l_w)
    xlim([-3 60])
    ylim([.77 1.025])
    set(f1,'Position',[10 10 l1_w h1_w]);
    
    f2=figure (num_fig+1);
    hold on
    grid on
    plot(time,A./S_tot,type_col_line,'LineWidth',l_w);
    xlim([-3 60])
    ylim([-0.025 0.5])
    set(f2,'Position',[10 10 l1_w h1_w]);
end


% load and plot the data

% load the volume data
name_data_3b='data_muzzey_09/fig_3b.txt';
str_data_3b=importdata(name_data_3b);
data_3b=str_data_3b.data;
data_3b_time=data_3b(:,1);

% we don't plot the data after the volume recovering ...
% in our model we don't consider any growth mechanism, so
% when the volume returns to its prestimulus value (V=1)
% then the followings data points are eqaul to 1

data_3b_time_v=data_3b(4:end,1);
data_3b_exp=ones(31,1);
data_3b_exp(1:15)=1+data_3b(4:18,6);
data_3b_std(1:15)=data_3b(4:18,7);

data_v_time=data_3b_time_v;
data_v_exp=data_3b_exp;
data_v_std=data_3b_std;

% load the Hog data
name_data_3a='data_muzzey_09/fig_3a.txt';
str_data_3a=importdata(name_data_3a);
data_3a=str_data_3a.data;
t_0=1;

data_3a_time=data_3a(t_0:end,1);
data_hog_exp=data_3a(t_0:end,6);
data_hog_std=data_3a(t_0:end,7);
  
f1=figure (num_fig);
hold on
grid on
errorbar(data_v_time(1:15),data_v_exp(1:15),data_v_std(1:15),'dr','MarkerSize',1,'LineWidth',0.1);
ylabel('Relative volume')
xlabel('time [min]')

f2=figure (num_fig+1);
hold on
grid on
ts=2;
errorbar(data_3a_time(1:ts:end),data_hog_exp(1:ts:end),data_hog_std(1:ts:end),...
    'dr','MarkerSize',1,'LineWidth',0.1);
ylabel('Hog1 activity')
xlabel('time [min]')

