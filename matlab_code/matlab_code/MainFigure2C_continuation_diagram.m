% The script generates the continuation diagram shown in Fig. 2C 
% for the two-state system - case 3  

clc
clear all


% nominal set
k=1;gamma=1;alpha=1;beta=1;k_u=1;
Kx=0.6; nx=4;


options=[];

%figure setting
line_w=1;
w_1=300;
w_2=300;

Ky=0.01:0.01:0.6;


ZZ12=[];
ZZ22=[];
ZZ32=[];

ZZ14=[];
ZZ24=[];
ZZ34=[];

ZZ18=[];
ZZ28=[];
ZZ38=[];

for idx_Ky=1:length(Ky)
    
    my=k*Ky(idx_Ky)-gamma/2;
    mx=alpha*Kx-beta/2;
    
    % disturbance/input of amplitude a_du (du1,du2,du3) applied at time tu
    % (tu1, tu2, tu3)
    
    tu1=0;
    du1=0;
    
    tu2=0;
    du2=0.1;
    
    tu3=0;
    du3=-0.1;
    
    % simulation time
    tspan=-20:0.1:100;
     
    % initial conditions
    y10= Ky(idx_Ky);
    y20= Kx;
    y0=[y10 y20]; 
    
    ny=2;
    
    [T12,Z12]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu1,du1,k_u,options);
    [T22,Z22]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu2,du2,k_u,options);
    [T32,Z32]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu3,du3,k_u,options);
    
    ZZ12=[ZZ12;Z12(end,:)];
    ZZ22=[ZZ22;Z22(end,:)];
    ZZ32=[ZZ32;Z32(end,:)];
    
    ny=4;
    
    [T14,Z14]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu1,du1,k_u,options);
    [T24,Z24]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu2,du2,k_u,options);
    [T34,Z34]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu3,du3,k_u,options);
    
    ZZ14=[ZZ14;Z14(end,:)];
    ZZ24=[ZZ24;Z24(end,:)];
    ZZ34=[ZZ34;Z34(end,:)];
    
    ny=8;
    
    [T18,Z18]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu1,du1,k_u,options);
    [T28,Z28]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu2,du2,k_u,options);
    [T38,Z38]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky(idx_Ky),tu3,du3,k_u,options);
    
    ZZ18=[ZZ18;Z18(end,:)];
    ZZ28=[ZZ28;Z28(end,:)];
    ZZ38=[ZZ38;Z38(end,:)];
    
end

fig=figure;
hold on
grid on
plot(Ky,ZZ12(:,1), 'k' ,Ky,ZZ22(:,1),'b--',Ky,ZZ32(:,1),'b-.','LineWidth',line_w);
plot(Ky,ZZ24(:,1),'r--',Ky,ZZ34(:,1),'r-.','LineWidth',line_w);
plot(Ky,ZZ28(:,1),'g--',Ky,ZZ38(:,1),'g-.','LineWidth',line_w);
xlabel('K_y')
ylabel('y')
axis([0.01 0.4 0 0.4])
grid on
legend('u_d=0','n_y=2, u_d^{+}','n_y=2, u_d^{-}','n_y=4, u_d^{+}','n_y=4, u_d^{-}','n_y=8, u_d^{+}','n_y=8, u_d^{-}')
set(fig,'Position',[10 10 w_1 w_2]);
set(gca,'FontSize',12)

