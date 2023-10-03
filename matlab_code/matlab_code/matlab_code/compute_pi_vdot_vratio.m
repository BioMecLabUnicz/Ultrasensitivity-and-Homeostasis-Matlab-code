function [PIi, PIe, PIt,  Vdot, Vratio]=compute_pi_vdot_vratio(T,x)
% compute the pressures, the volume and its variation


% parameters for the glycerol and biophysical modules
global global_V_b  global_PIi_0 global_PIe_0 global_glyc_0;
global global_kp1  global_V_0 global_V_PIt0;
global global_t0 global_PIe_d;



kp1=global_kp1; % bounds [5.2e-3, 160]
V_0=global_V_0;
V_PIt0=global_V_PIt0; % bounds [ 0.5, 0.99]
V_b=global_V_b; %[0.31 0.46]
PIi_0=global_PIi_0; %[0.6 0.7]
PIe_0=global_PIe_0; % [0.24 0.25]
glyc_0=global_glyc_0; %[1.1e-4 5e-4]
t0=global_t0;
PIe_d=global_PIe_d;
PIt_0=PIi_0-PIe_0;
n=PIi_0*(V_0-V_b)-glyc_0;

%outputs
V=x(:,9);
glyc=x(:,10);

% internal pressure
PIi=(n+glyc)./(V-V_b);

% external pressure
PIe=zeros(length(T),1);
for t=1:length(T)
    
    if length(t0)==1
        if T(t)<t0
            PIe(t)=PIe_0;
        else
            PIe(t)=PIe_0+PIe_d;
        end
        
    else
        if T(t)<t0(1)
            PIe(t)=PIe_0;
        elseif T(t)>=t0(end)
            PIe(t)=PIe_0+PIe_d(end);
        else
            k=1;
            t_found=false;
            while not(t_found)
                if T(t)>=t0(k) && T(t)<t0(k+1)
                    t_found=true;
                    PIe(t)=PIe_0+PIe_d(k);
                end
                k=k+1;
            end
        end
    end
end


% turgor pressure
PIt=zeros(length(T),1);
for t=1:length(T)
    if V(t)>V_PIt0
        PIt(t)=PIt_0.*(V(t)-V_PIt0)/(V_0-V_PIt0);
    else
        PIt(t)=0;
    end
end

% ode for Volume
Vdot=kp1.*(PIi-PIe-PIt);

Vratio=Vdot./(V-V_b);


