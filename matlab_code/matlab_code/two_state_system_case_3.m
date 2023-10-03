% ODE model of the two-state system for case 3
function dy=two_state_system_case_3(t,y,my,k,gamma,nx,Kx,mx,alpha,beta,ny,Ky,tu,du,k_u,options)

dy=zeros(2,1);

%noise on the input - i.e. perturbation on the system 
if isempty(tu)    
    % no distrurb
    u_d=0;
else
    % disturb at time t0
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
            kk=1;
            t_found=false;
            while not(t_found)
                if t>=tu(kk) && t<tu(kk+1)
                    t_found=true;
                    u_d=du(kk);
                end
                kk=kk+1;
            end
        end
    end
end



dy(1)=my-k*y(1)+gamma*(y(2)^nx/(y(2)^nx+Kx^nx))-k_u*u_d;
dy(2)=mx-alpha*y(2)+beta*(Ky^ny/(y(1)^ny+Ky^ny));


return