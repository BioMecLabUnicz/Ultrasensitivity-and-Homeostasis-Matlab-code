% The script generates the time responses (shown in Fig. 2B) of the
% two-state system for case 2 to the applied step disturbance, 
% with Ky = 0.4 and = 0.1 for different values of ny 
% (ny=2 (blue curves), ny=4 (red) and ny=8 (green)). 

clc
clear all
close all

% vector for Ky values
v_Ky=[0.4  0.1];
% vector for ny values
v_ny=[2 4 8];

% step disturbance of amplitude a_du   
a_ud=0.1;
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
ky=1;gamma=1;alpha=1;beta=1;mx=1; k_u=1;

for idx_Ky=1:length(v_Ky)
    Ky=v_Ky(idx_Ky); % Ky value
    my=(ky*alpha*Ky-gamma*(mx+beta/2))*(1/alpha);
    
    % initial conditions
    y10= Ky;
    y20= (ky*Ky-my)/gamma;
    y0=[y10 y20]; 
    
    % step disturbance of amplitude a_du applied at time tu
    tu=[10 20 30];
    du=[a_ud 0 -a_ud];
    
    if idx_plot_var==1
        var_norm=y10; % plot y, nomalize if tag_y_norm=1
    elseif idx_plot_var==2
        var_norm=y20; % plot x, nomalize if tag_y_norm=1
    end
    if tag_subplot
        subplot(length(v_Ky),1,idx_Ky)
    else
        figure
        hold on
        grid on
    end
    
    
    for idx_ny=1:length(v_ny)
        ny=v_ny(idx_ny); % ny value
        
        % compute time response of the two-state system for case 2
        [T0,Z0_ny]=ode15s(@two_state_system_case_2,tspan,y0,[],my,ky,gamma,Ky,alpha,beta,mx,ny,tu,du,k_u,[]);

        if ny==2
            type_line_col='-b';
        elseif ny==4
            type_line_col='-r';
        elseif ny==8
            type_line_col='-g';
        else
            type_line_col='-k';
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
                ymin=0.34;
                ymax=0.46;
            elseif Ky==0.1
                ymin=0.06;
                ymax=0.14;
            end
        end
        ylim([ymin ymax])
        ylabel('y')
        if idx_ny==length(v_Ky)
            xlabel('time')
        end
        str_title=['K_y=',num2str(Ky)];
        title(str_title)
        set(gca,'FontSize',12)
    end    
end
set(f,'Position',[10 10 w_1 w_2]);
