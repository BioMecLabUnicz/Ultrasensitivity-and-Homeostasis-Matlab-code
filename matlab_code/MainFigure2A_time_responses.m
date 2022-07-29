% The script generates the time responses (shown in Fig. 2A) of the 
% two-state system for case 1 to the applied step disturbance, 
% with Kx = 0.4 and = 0.1 for different values of nx 
% (nx=2 (blue curves), nx=4 (red) and nx=8 (green)). 


clc
clear all
close all

% vector for Kx values
v_Kx=[ 0.4 0.1];
% vector for nx values
v_nx=[2 4 8];

% step disturbance of amplitude a_du   
a_ud=.1;
tag_y_norm=0; % if one, normalize the output, if zero no
idx_plot_var=1; %plot y
%idx_plot_var=2; %plot x

% simulation time
tspan=0:0.01:40;

% figure setting
l_w=1;
w_1=300;
w_2=300;
tag_subplot=1; %if one use subplot
if tag_subplot
    f=figure;
end

% define parameters
ky=1;gamma=1;alpha=1;beta=1;mx=1;k_u=1;

for idx_Kx=1:length(v_Kx)
    Kx=v_Kx(idx_Kx);  % Kx value
   
    my=(-alpha*ky*Kx-beta*gamma/2+ky*mx)/beta;

    % initial conditions

    y10=(my+gamma/2)/ky;
    y20=Kx;
    p1=[y10 y20];
    y0=p1;
    
    % step disturbance of amplitude a_du applied at time tu 
    tu=[10 20 30];
    du=[a_ud 0 -a_ud];
    
    if idx_plot_var==1
        var_norm=y10; % plot y, nomalize if tag_y_norm=1
    elseif idx_plot_var==2
        var_norm=y20; % plot x, nomalize if tag_y_norm=1
    end
    
    if tag_subplot
        subplot(length(v_Kx),1,idx_Kx)
    else
        f=figure;
        hold on
        grid on
    end
    for idx_nx=1:length(v_nx)
        nx=v_nx(idx_nx); % nx value
        
        % compute time response of the two-state system for case 1
        [T0,Z0_nx]=ode15s(@two_state_system_case_1,tspan,y0,[],my,ky,gamma,Kx,alpha,beta,mx,nx,tu,du,k_u,[]);

        if nx==2
            type_line_col='-b';
        elseif nx==4
            type_line_col='-r';
        elseif nx==8
            type_line_col='-g';
        else
            type_line_col='-k';
        end
        
        hold on 
        grid on
        if tag_y_norm
            plot(tspan,Z0_nx(:,idx_plot_var)/var_norm,type_line_col,'linewidth',l_w);
            ymin=0.75;
            ymax=1.25;
        else
            plot(tspan,Z0_nx(:,idx_plot_var),type_line_col,'linewidth',l_w); 
            if Kx==0.4
                ymin=0.54;
                ymax=0.66;
            elseif Kx==0.1
                ymin=0.86;
                ymax=0.94;
            end
        end
        ylim([ymin ymax])
        ylabel('y')
        if idx_nx==length(v_Kx)
            xlabel('time')
        end
         str_title=['K_x=',num2str(Kx)];
        title(str_title)
        set(gca,'FontSize',12)
    end
end
set(f,'Position',[10 10 w_1 w_2]);