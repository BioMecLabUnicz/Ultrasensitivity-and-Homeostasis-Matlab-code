% The script generates the results shown in Figure 3A  (or Figure S3A and Figure S3B) 
% by setting Kx=0.4 or 0.1 and by varying the parameters between 10% or (50%) of their nominal values. 
% You set the number of simulations (num_sim) by the input comand (in the
% ms the results are obtained by setting num_sim=1000)


clc
clear all
close all

num_sim=input('number of simulations: ');


% varying the parameters between +-10% of their nominal values
per=0.1;
% varying the parameters between +-50% of their nominal values
%per=0.5;

% set the value for Kx, 
Kx=0.4; % case 1 with K_x=0.4
%Kx=0.1; % case 1 with K_x=0.1

% nx values
v_nx=[2 4 8];

options=[];
% nominal set
k=1;gamma=1;alpha=1;beta=1;mx=1;k_u=1;
my=(-alpha*k*Kx-beta*gamma/2+k*mx)/beta;

% simulation time
tspan=-10:0.1:40;

% initial conditions
y01=(my+gamma/2)/k;
y02=Kx;
y0=[y01 y02]; % initial equil.

% disturbance/input of amplitude a_du applied at time tu
tu=[10 20 30];
perc_dist=0.3;

a_ud0=y01*perc_dist*k;
du0=[a_ud0 0 -a_ud0];


% parameters for generating figures
array_no={'','',''};
num_font=10;
array_n={'n_x=2','n_x=4','n_x=8'};
l_w=600;
h_w=500;
l_w1=550;
h_w1=450;

y_low=0.69;
y_high=1.31;

box_low=-0.01;
box_up=0.31;

line_w=1;

fa=figure;
hold on
grid on

length_v_nx=length(v_nx);

err_n_pos=NaN(1,length_v_nx);
err_n_neg=NaN(1,length_v_nx);

err_19=NaN(length_v_nx,num_sim);
err_39=NaN(length_v_nx,num_sim);


for idx_nx=1:length_v_nx
    Z=[];
    Z_scale=[];
    ZZ=[];
    ZZ_scale=[];
    
    if idx_nx==1       
        PP=[]; % matrix used for saving at each simulation the parameters which are 
                  % varied randomly between +/-10% (or 50%) of their nominal values
    else
        eval (['load sys_case1','_Kx_',num2str(Kx*100),'_distur_',num2str(perc_dist*100),'_perc100_par_u1']);
    end
     
    nx=v_nx(idx_nx);
    
    for idx_sim=1:num_sim
            
        if idx_nx==1
            
            %in general, you can generate N random numbers in the interval (a,b)
            % with the formula r = a + (b-a).*rand(N,1).
        
            k1 = k-per*k+(k+per*k-(k-per*k))*rand(1,1);
            Kx1=Kx-per*Kx+(Kx+per*Kx-(Kx-per*Kx))*rand(1,1);
            gamma1 = gamma-per*gamma+(gamma+per*gamma-(gamma-per*gamma))*rand(1,1);
            mx1 = mx-per*mx+(mx+per*mx-(mx-per*mx))*rand(1,1);
            alpha1 =  alpha-per*alpha+(alpha+per*alpha-(alpha-per*alpha))*rand(1,1);
            beta1 =  beta-per*beta+(beta+per*beta-(beta-per*beta))*rand(1,1);
            
            my0=my;          
            my1 = my0-per*my0 + (my0+per*my0-(my0-per*my0))*rand(1,1);
            
            
           
            [~,Zn1]=ode15s(@two_state_system_case_1,tspan,y0,[],my1,k1,gamma1,Kx1,alpha1,beta1,mx1,v_nx(1),[],[],k_u,options);               
            y_eq1=Zn1(end,1);
            x_eq1=Zn1(end,2);
            [~,Zn2]=ode15s(@two_state_system_case_1,tspan,y0,[],my1,k1,gamma1,Kx1,alpha1,beta1,mx1,v_nx(2),[],[],k_u,options);   
            y_eq2=Zn2(end,1);
            x_eq2=Zn2(end,2);
            [~,Zn3]=ode15s(@two_state_system_case_1,tspan,y0,[],my1,k1,gamma1,Kx1,alpha1,beta1,mx1,v_nx(3),[],[],k_u,options);
            y_eq3=Zn3(end,1);
            x_eq3=Zn3(end,2);
            if  (y_eq1 >0) &&  (y_eq2 >0) &&  (y_eq3 >0) && (x_eq1 >0) &&  (x_eq2 >0) &&  (x_eq3 >0)
                min_yeq=min([y_eq1 y_eq2 y_eq3]);
                
                a_ud1=min_yeq*perc_dist*k1;
            else
                a_ud1=0;
            end
                
            P=[my1; k1; gamma1; Kx1; mx1; alpha1; beta1;a_ud1];



            PP=cat(2,PP,P);
            
            
            
        else
            
            my1 = PP(1,idx_sim);
            k1 = PP(2,idx_sim);
            gamma1 = PP(3,idx_sim);
            Kx1=PP(4,idx_sim);
            mx1 = PP(5,idx_sim);
            alpha1 =  PP(6,idx_sim);
            beta1 =  PP(7,idx_sim); 
            a_ud1 = PP(8,idx_sim);
            
        end
        
        % using the varied parameter set
     
        du=[a_ud1 0 -a_ud1];
        [T1,Z1]=ode15s(@two_state_system_case_1,tspan,y0,[],my1,k1,gamma1,Kx1,alpha1,beta1,mx1,nx,tu,du,k_u,options);
        time9=find(tspan==9);
        time19=find(tspan==19);
        time39=find(tspan==39);
        y_eq=Z1(time9,1);
        x_eq=Z1(time9,2);
        y_eq1=Z1(time19,1);        
        y_eq2=Z1(time39,1);

        x_eq1=Z1(time19,2);        
        x_eq2=Z1(time39,2);

    
                
        if y_eq<0|| x_eq<0 || y_eq1<0|| x_eq1<0 || y_eq2<0|| x_eq2<0 || a_ud1==0
           
            Z1=NaN(length(T1),2); 
        end
        Z=cat(2,Z,Z1(:,1));
        ZZ=cat(2,ZZ,Z1(:,2));
        Z_scale=cat(2,Z_scale,Z1(:,1)/y_eq);
        ZZ_scale=cat(2,ZZ_scale,Z1(:,2)/x_eq);
    
    end
    
    if idx_nx==1
        eval (['save sys_case1','_Kx_',num2str(Kx*100),'_distur_',num2str(perc_dist*100),'_perc100_par_u1   PP']);
    end 
    eval (['save sys_case1','_Kx_',num2str(Kx*100),'_nx_',num2str(nx),'_distur_',num2str(perc_dist*100),'_perc100_par_u1  Z Z_scale ZZ ZZ_scale PP']);
        
    
    time9=find(tspan==9);
    time19=find(tspan==19);
    time39=find(tspan==39);
     
    % compute the error as defined in the ms
    err_19(idx_nx,:)=abs(Z_scale(time9,:)-Z_scale(time19,:));
    err_39(idx_nx,:)=abs(Z_scale(time9,:)-Z_scale(time39,:));
  
    
    % using the nominal parameter set
    [T0,Z0_nx]=ode15s(@two_state_system_case_1,tspan,y0,[],my,k,gamma,Kx,alpha,beta,mx,nx,tu,du0,k_u,[]);

    if idx_nx==1
        type_lin_col='-b';
    elseif idx_nx==2
        type_lin_col='-r';
    elseif idx_nx==3
        type_lin_col='-g';
    else
        type_lin_col='-c';
    end
    
    % compute the error for the nominal case 
    err_n_pos(idx_nx)=abs(Z0_nx(time9,1)/y01-Z0_nx(time19,1)/y01);
    err_n_neg(idx_nx)=abs(Z0_nx(time9,1)/y01-Z0_nx(time39,1)/y01);

   
%     subplot(length(v_nx),1,idx_nx)

    if idx_nx==1
        subplot(6,2,[1 3])
    elseif idx_nx==2
        subplot(6,2,[5 7])
     elseif idx_nx==3
        subplot(6,2,[9 11])
    end
    
    hold on
    grid on
    plot(tspan,Z_scale(:,:),'-k','linewidth',line_w)
    plot(tspan,Z0_nx(:,1)/y01,type_lin_col,'linewidth',line_w)

    xlim([0 40])
    ylabel('y_{norm}')
    ylim([y_low y_high])
    if idx_nx==3
        xlabel('time')
    end
    str_title=['n_x=',num2str(nx)];
    title(str_title)

    set(gca,'FontSize',num_font)

end


 
subplot(6,2,[2 4 6]) 
hold on
grid on

boxplot(err_19','labels',array_n)
plot(1.4,err_n_pos(1),'*b','markersize',4)
plot(2.4,err_n_pos(2),'*r','markersize',4)
plot(3.4,err_n_pos(3),'*g','markersize',4)
ylabel('e_n^+')
ylim([box_low box_up])
set(gca,'FontSize',num_font)
set(gca,'TickLabelInterpreter','tex')

subplot(6,2,[8 10 12]) 
hold on
grid on
boxplot(err_39','labels',array_n)
plot(1.4,err_n_neg(1),'*b','markersize',4)
plot(2.4,err_n_neg(2),'*r','markersize',4)
plot(3.4,err_n_neg(3),'*g','markersize',4)
ylabel('e_n^-')
ylim([box_low box_up])

set(gca,'FontSize',num_font) 
set(gca,'TickLabelInterpreter','tex')
%set(fa,'Position',[10 10 l_w1 h_w1 ]);
