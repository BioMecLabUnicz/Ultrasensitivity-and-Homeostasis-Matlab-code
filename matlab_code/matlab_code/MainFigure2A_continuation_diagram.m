% The script generates the continuation diagram shown in Fig. 2A 
% for the two-state system - case 1.  
clc
clear all

k_y=1;gamma=1;alpha=1;beta=1;mx=1;k_u=1;

options=[];

%figure setting
line_w=1;
w_1=300;
w_2=300;

Kx=0.01:0.01:0.4;


ZZ12=[];
ZZ22=[];
ZZ32=[];

ZZ14=[];
ZZ24=[];
ZZ34=[];

ZZ18=[];
ZZ28=[];
ZZ38=[];

for idx_Kx=1:length(Kx)
    my=(-alpha*k_y*Kx(idx_Kx)-beta*gamma/2+k_y*mx)/beta;
    
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
    y10=(my+gamma/2)/k_y;
    y20=Kx(idx_Kx);
    p1=[y10 y20];
    y0=p1;
    
    
    nx=2;
    
    [T12,Z12]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu1,du1,k_u,options);
    [T22,Z22]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu2,du2,k_u,options);
    [T32,Z32]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu3,du3,k_u,options);
    
    ZZ12=[ZZ12;Z12(end,:)];
    ZZ22=[ZZ22;Z22(end,:)];
    ZZ32=[ZZ32;Z32(end,:)];
    
    nx=4;
    
    [T14,Z14]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu1,du1,k_u,options);
    [T24,Z24]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu2,du2,k_u,options);
    [T34,Z34]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu3,du3,k_u,options);
    
    ZZ14=[ZZ14;Z14(end,:)];
    ZZ24=[ZZ24;Z24(end,:)];
    ZZ34=[ZZ34;Z34(end,:)];
    
    
    nx=8;
    
    [T18,Z18]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu1,du1,k_u,options);
    [T28,Z28]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu2,du2,k_u,options);
    [T38,Z38]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k_y,gamma,Kx(idx_Kx),alpha,beta,mx,nx,tu3,du3,k_u,options);
    
    ZZ18=[ZZ18;Z18(end,:)];
    ZZ28=[ZZ28;Z28(end,:)];
    ZZ38=[ZZ38;Z38(end,:)];
    
end


fig=figure;
hold on
grid on
plot(Kx,ZZ12(:,1), 'k' ,Kx,ZZ22(:,1),'b--',Kx,ZZ32(:,1),'b-.','LineWidth',line_w);
plot( Kx,ZZ24(:,1),'r--',Kx,ZZ34(:,1),'r-.','LineWidth',line_w);
plot( Kx,ZZ28(:,1),'g--',Kx,ZZ38(:,1),'g-.','LineWidth',line_w);
legend('u_d=0','n_x=2, u_d^{+}','n_x=2, u_d^{-}','n_x=4, u_d^{+}','n_x=4, u_d^{-}','n_x=8, u_d^{+}','n_x=8, u_d^{-}')
xlabel('K_x')
ylabel('y')
axis([0.01 0.4 0.55 1])
set(fig,'Position',[10 10 w_1 w_2]);
set(gca,'FontSize',12)
