% The script generates the isoclines shown in Fig. 2C for the two-state system - case 3 -
% by setting K_y (=0.4 and 0.1 as shown in Fig. 2C of the ms).

clear all
clc

% choose K_y value
K_y=input('set the value for K_y in [0.1 0.5] '); 
%K_y=.1; 
%K_y=.4; 

fig=figure;
hold on
grid on

% figure setting
w_1=300;
w_2=300;
line_w=1;

% y vector
y_ss=0:0.01:1;

% parameters
n_y_v=[2 4 8];

K_x=0.6;
n_x=4;
u=0.1;
alpha=1;
beta=1;
gamma=1; 
k=1; 

for idx_n_y=1:length(n_y_v)
        
    % n_y value
    n_y=n_y_v(idx_n_y);    
    
    m_y=k*K_y-gamma/2;
    m_x=alpha*K_x- beta/2;

    f_y_ss=K_y^n_y./(K_y^n_y+y_ss.^n_y);
    x_ss_nonlin=(m_x+beta.*f_y_ss)./alpha;
    x_ss_contr=x_ss_nonlin.^n_x./(K_x^n_x + x_ss_nonlin.^n_x);
        
    r0=k*y_ss-m_y;
    r2=k*y_ss-m_y+u;
    r1=k*y_ss-m_y-u;
    
    if idx_n_y==1
        tag_line_col='-b';
    elseif idx_n_y==2
        tag_line_col='-r';
    elseif idx_n_y==3
        tag_line_col='-g';
    else
        tag_line_col='-c';
    end
    
    % plot h (with different n_y values) vs y
    plot(y_ss,gamma*x_ss_contr,tag_line_col,'linewidth',line_w);
    
end

% plot r for different values of the disturbance amplitude
plot(y_ss,r0,'-k','linewidth',line_w)
plot(y_ss,r1,'-.k','linewidth',line_w)
plot(y_ss,r2,'--k','linewidth',line_w)
ylim([0 1])
xlim([0 .65])
legend('h, n_y=2','h, n_y=4','h, n_y=8','r with u_d=0','r with u_d=0.1','r with u_d=-0.1')
str_title=['K_y=',num2str(K_y)];
title(str_title)
ylabel('h, r')
xlabel('y')
set(fig,'Position',[10 10 w_1 w_2]);
set(gca,'FontSize',12)
