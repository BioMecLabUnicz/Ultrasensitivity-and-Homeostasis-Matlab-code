% The script generates the time responses (shown in Fig. S1) of the
% two-state system for case 3 to the applied step disturbance, 
% with Ky = 0.6 or 0.4 and Kx=0.6 or 0.4  for different values of nx 
% (nx=2 (blue curves), nx=4 (red) and nx=8 (green)) and two different
% values for ny (ny=2, solid lines; ny=4, dashed lines)

clc
clear all
close all


%choose Kx and Ky values
Kx=.6;
%Kx=.4;

Ky=.6;
%Ky=.4;

% vector for ny values
v_ny=[2 4 ];
% vector for nx values,
v_nx=[2 4 8];

% step disturbance of amplitude a_du   
a_ud=.1;

idx_plot_var=1; %plot y
%idx_plot_var=2; %plot x
tag_y_norm=0; % if one, normalize the output, if zero no

% figure setting
l_w=2;
w_1=450;
w_2=300;

tag_subplot=0; %if one use subplot
fig=figure;
title_plot=1;
hold on
grid on

% define parameters
k=1;gamma=1;alpha=1;beta=1;k_u=1;

% simulation time
tspan=0:0.01:40;

for idx_ny=1:length(v_ny)
    
    
    ny=v_ny(idx_ny);
    my=k*Ky-gamma/2;
    mx=alpha*Kx- beta/2;
    
    % intial conditions
    y10= Ky;
    y20= Kx;
    y0=[y10 y20]; 
    
    % step disturbance of amplitude a_du applied at time tu 
    tu=[10 20 30];
    du=[a_ud 0 -a_ud];

    if idx_plot_var==1
        var_norm=y10; % plot y, nomalize if tag_y_norm=1
    elseif idx_plot_var==2
        var_norm=y20;  % plot x, nomalize if tag_y_norm=1
    end
    
    if tag_subplot
        subplot(length(v_ny),1,idx_ny)
    
    end
    
    for idx_nx=1:length(v_nx)
    
        % nx valaue
        nx=v_nx(idx_nx);
        % compute time response of the two-state system for case 3
        [T0,Z0_ny]=ode15s(@two_state_system_case_3,tspan,y0,[],my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky,tu,du,k_u,[]);
        if ny==2
            if nx==2
                type_line_col='-b';
            elseif nx==4
                type_line_col='-r';
            elseif nx==8
                type_line_col='-g';
            else
                type_line_col='-k';
            end
        elseif ny==4

            if nx==2
                type_line_col='-.b';
            elseif nx==4
                type_line_col='-.r';
            elseif nx==8
                type_line_col='-.g';
            else
                type_line_col='-.k';
            end
        end
        hold on
        grid on
        
        if tag_y_norm
            plot(tspan,Z0_ny(:,idx_plot_var)/var_norm,type_line_col,'linewidth',l_w);
            ymin=0.8;
            ymax=1.2;
        else
            plot(tspan,Z0_ny(:,idx_plot_var),type_line_col,'linewidth',l_w);
            if Ky==0.4
                ymin=0.345;
                ymax=0.455;
            elseif Ky==0.2
                ymin=0.16;
                ymax=0.24;
            elseif Ky==0.1
                ymin=0.06;
                ymax=0.14;
                
            elseif Ky==0.6
                ymin=0.53;
                ymax=0.67;
            end
        end
        ylim([ymin ymax])
        ylabel('y')
        xlabel('time')
        if title_plot
            str_title=['K_x=',num2str(Kx), '; K_y=',num2str(Ky),];
            title(str_title)
        end
        set(gca,'FontSize',14)
    end       
end
legend('n_x=2, n_y=2','n_x=4, n_y=2','n_x=8, n_y=2','n_x=2, n_y=4','n_x=4, n_y=4','n_x=8, n_y=4')
set(fig,'Position',[10 10 w_1 w_2]);

