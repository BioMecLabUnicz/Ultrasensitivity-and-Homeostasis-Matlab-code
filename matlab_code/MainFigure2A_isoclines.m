% The script generates the isoclines shown in Fig. 2A for the two-state system - case 1 -
% by setting K_x (=0.4 and 0.1 as shown in Fig. 2A of the ms).

clc
clear all

% choose K_x value
K_x=input('set the value for K_x in [0.1 0.5] '); 

% y vector
y_ss=0:0.01:1;

% parameters
alpha=1;
beta=1;
gamma=1; 
k=1; 
m_x=1;    
m_y=k*(m_x-K_x*alpha)/beta-gamma/2; % given a value for K_x
y_eq=(m_y+gamma/2)/k;
u=0.1; 

fig=figure;
hold on
grid on

w_1=300;
w_2=300;

line_w=1;
n_x_v=[2 4 8];
   
for idx_n_x=1:length(n_x_v)
    
    % n_x value
    n_x=n_x_v(idx_n_x);
     
    x_ss_lin=(m_x-beta.*y_ss)./alpha;
    x_ss_contr_lin=gamma*x_ss_lin.^n_x./(K_x^n_x + x_ss_lin.^n_x);
    
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
    plot(y_ss,x_ss_contr_lin,tag_line_col,'linewidth',line_w);
    
end

% plot r for different values of the disturbance amplitude
plot(y_ss,r0,'-k','linewidth',line_w')
plot(y_ss,r2,'--k','linewidth',line_w')
plot(y_ss,r1,'-.k','linewidth',line_w)
ylim([0 1])
legend('h, n_x=2','h, n_x=4','h, n_8=8','r_0','r_{u^{+}}','r_{u^{-}}')
ylabel('h, r')
xlabel('y')
set(fig,'Position',[10 10 w_1 w_2]);
set(gca,'FontSize',12)
