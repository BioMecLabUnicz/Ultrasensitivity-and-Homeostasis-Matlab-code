% The script generates the isoclines shown in Fig. S1 for the two-state system - case 3 -
% with different f function values (K_x=0.6 or 0.4 with n_x=2, 4, 8) and g
% values (K_y =0.6 or 0.4, ny=2 or 4)

clear all
clc
close all

%choose K_x value
K_x=.6;
%K_x=.4;

% choose K_y and n_y values
K_y=.6;
% K_y=.4;

n_y=2;
%n_y=4;    


fig=figure;
hold on
grid on
title_subplot=1;

% figure setting
w_1=300;
w_2=300;
line_w=2;

% y vector
y_ss=0:0.01:1;

% parameters
n_x_v=[2 4 8];

u=0.1;
alpha=1;
beta=1;
gamma=1; 
k=1; 

for idx_n_x=1:length(n_x_v)
        
    % n_x value
    n_x=n_x_v(idx_n_x);
    m_y=k*K_y-gamma/2;
    m_x=alpha*K_x- beta/2;

    f_y_ss=K_y^n_y./(K_y^n_y+y_ss.^n_y);
    x_ss_nonlin=(m_x+beta.*f_y_ss)./alpha;
    x_ss_contr=x_ss_nonlin.^n_x./(K_x^n_x + x_ss_nonlin.^n_x);
        
    r0=k*y_ss-m_y;
    r2=k*y_ss-m_y+u;
    r1=k*y_ss-m_y-u;
    
    if idx_n_x==1
        tag_line_col='-b';
    elseif idx_n_x==2
        tag_line_col='-r';
    elseif idx_n_x==3
        tag_line_col='-g';
    else
        tag_line_col='-c';
    end
    
    % plot h (with different n_x values) vs y
    plot(y_ss,gamma*x_ss_contr,tag_line_col,'linewidth',line_w);
    
end

% plot r for different values of the disturbance amplitude
plot(y_ss,r0,'-k','linewidth',line_w)
plot(y_ss,r1,'-.k','linewidth',line_w)
plot(y_ss,r2,'--k','linewidth',line_w)
ylim([0 1])
xlim([0 0.8])
legend('h, n_x=2','h, n_x=4','h, n_x=8','r_0','r_{u_d^+}','r_{u_d^-}')
if title_subplot
    str_title=['K_x=',num2str(K_x), '; K_y=',num2str(K_y), ', n_y=',num2str(n_y)];
    title(str_title)
end

ylabel('isoclines h, r')
xlabel('y')
set(fig,'Position',[10 10 w_1 w_2]);
set(gca,'FontSize',14)
