% ODE model of the two-state system for case 2
function dy=two_state_system_case_2(t,y,my,ky,gamma,Ky,alpha,beta,mx,ny,tu,du,k_u,options)
         
dy=zeros(2,1);

%noise on the input - i.e. perturbation on the system 
if isempty(tu)    
    % no distrurb
    u_d=0;
else
    % disturb at time tu(i)
    if length(tu)==1
        if t<tu
            u_d=0;
        else
            u_d=du;
        end
    else
        if t<tu(1)
             u_d=0;
        elseif t>=tu(end)
            u_d=du(end);
        else
            k=1;
            t_found=false;
            while not(t_found)
                if t>=tu(k) && t<tu(k+1)
                    t_found=true;
                    u_d=du(k);
                end
                k=k+1;
            end
        end
    end
end

dy(1)=-ky*y(1)+gamma*y(2)+my-k_u*u_d;
dy(2)=-alpha*y(2)+beta*(Ky^ny/(y(1)^ny+Ky^ny))+mx;


return