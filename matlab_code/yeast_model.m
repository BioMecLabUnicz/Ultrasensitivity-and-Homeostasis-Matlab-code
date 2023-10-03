% ODE model for yeast osmosensing system

function xdot=yeast_model(t,x,flag)
                        

% parameters for the SLN1 sensor
global ka ks krs km krm khs khm 

global krc kc krt kt khc

% parameters for the HOG MAPK
global a1 a2  k1 k2 d1 d2

% parameters for the glycerol and biophysical modules
global global_V_b global_V_e global_PIi_0 global_PIe_0 global_glyc_0;
global global_kp1 global_kfps1 global_k_hog global_V_0 global_V_PIt0;
global global_t0 global_PIe_d;

global  n_hog_fps1 n_hog_glyc K_hog_glyc  K_hog_fps1  


kp1=global_kp1; % bounds [5.2e-3, 160]
kp2=global_kfps1;  % bounds (0, inf)
k_hog=global_k_hog;   % bounds (0, inf)
V_0=global_V_0;
V_PIt0=global_V_PIt0; % bounds [ 0.5, 0.99]
V_b=global_V_b; %[0.31 0.46]
V_e=global_V_e; %[0.5e3 5e3]
PIi_0=global_PIi_0; %[0.6 0.7]
PIe_0=global_PIe_0; % [0.24 0.25]
glyc_0=global_glyc_0; %[1.1e-4 5e-4]
t0=global_t0;
PIe_d=global_PIe_d;
PIt_0=PIi_0-PIe_0;
n=PIi_0*(V_0-V_b)-glyc_0;


% state variables

HKp=x(1);
RCp=x(2);
HPTp=x(3);
R1p=x(4);
R2p=x(5);

A=x(6);       % Free Active Protein Concentration
C1=x(7);      % Concentration of Inactive Enzyme-Substrate Complex, C1
C2=x(8);      % Concentration of Active Enzyme-Substrate Complex, C2


V=x(9);
glyc=x(10);
glyc_ex=x(11);


HK_tot=x(12);

RC_tot=x(13);

HPT_tot=x(14);

Rs_tot=x(15);

Rm_tot=x(16);

S_tot=x(17);

E2_tot=x(18);




% conservation

HK=HK_tot-HKp;

RC=RC_tot-RCp;

HPT=HPT_tot-HPTp;

R1=Rs_tot-R1p-C1;

R2=Rm_tot-R2p;

I=S_tot-A-C1-C2;
E2=E2_tot-C2;

%compute the PIt
if V>V_PIt0
    PIt=PIt_0*(V-V_PIt0)/(V_0-V_PIt0);
else
    PIt=0;
end

%compute the PIe
if length(t0)==1
    if t<t0
        PIe=PIe_0;
    else
        PIe=PIe_0+PIe_d;
    end       
else
    if t<t0(1)
        PIe=PIe_0;
    elseif t>=t0(end)
        PIe=PIe_0+PIe_d(end);
    else
        k=1;
        t_found=false;
        while not(t_found)
            if t>=t0(k) && t<t0(k+1)
                t_found=true;
                PIe=PIe_0+PIe_d(k);
            end
            k=k+1;
        end
    end
end

%compute the PIi
PIi=(n+glyc)/(V-V_b);


% ode for Volume
Vdot=kp1*(PIi-PIe-PIt);

Vratio=Vdot/(V-V_b);

c_ratio=1;
HK_tot_dot=-c_ratio*Vratio*HK_tot;

RC_tot_dot=-c_ratio*Vratio*RC_tot;

HPT_tot_dot=-c_ratio*Vratio*HPT_tot;

Rs_tot_dot=-c_ratio*Vratio*Rs_tot;

Rm_tot_dot=-c_ratio*Vratio*Rm_tot;

S_tot_dot=-c_ratio*Vratio*S_tot;

E2_tot_dot=-c_ratio*Vratio*E2_tot;


% Fps1 channel

Af=A;
Ah=A;
u_Fps1=kp2*(K_hog_fps1^n_hog_fps1)/(Af^n_hog_fps1+K_hog_fps1^n_hog_fps1);
u_Diff=u_Fps1*(glyc/(V-V_b)-glyc_ex/V_e);


% ode for external glycerol
glyc_exdot=u_Diff;

%glycerol production by hog
u_hog=k_hog*(Ah^n_hog_glyc)/(Ah^n_hog_glyc+K_hog_glyc^n_hog_glyc);

% ode for internal glycerol
glycdot=u_hog-u_Diff -(c_ratio*Vratio)*glyc;


% ODEs for phos  
dot_HKp  = ka*PIt*HK + krc*RCp*HK - kc*RC*HKp;

dot_RCp  = kc*RC*HKp + krt*HPTp*RC - krc*RCp*HK - kt*HPT*RCp - khc*RCp;
dot_HPTp = kt*HPT*RCp + krs*R1p*HPT + krm*R2p*HPT - krt*HPTp*RC - ks*R1*HPTp - km*R2*HPTp;
dot_R1p  = ks*R1*HPTp - krs*R1p*HPT - khs*R1p;
dot_R2p  = km*R2*HPTp - krm*R2p*HPT - khm*R2p;

% ODEs for mapk/hog
dA = k1*C1 - a2*A*E2 + d2*C2;       % A dot
dC1 = a1*I*R1 - d1*C1 - k1*C1;       % C1 dot
dC2= a2*E2*A - k2*C2 - d2*C2;       % C2 dot
    


xdot=[dot_HKp;dot_RCp;dot_HPTp;dot_R1p;dot_R2p;dA; dC1; dC2;Vdot;glycdot;...
    glyc_exdot;HK_tot_dot;RC_tot_dot;HPT_tot_dot;Rs_tot_dot;Rm_tot_dot;...
    S_tot_dot;E2_tot_dot];
