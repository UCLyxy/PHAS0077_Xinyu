function r=calculate_residual(p,u)
% calculate_residual: return residual(rhs) for PDE and auxiliary equations
%
%  USAGE:
%   r=calculate_residual(p,u)
%
%  INPUTS:
%   *p: object oriented pde object
%   *u: solution vector of pde system(u=[u1;u2;par])
%
%  RETURNS:
%   *r: residual for PDE and auxiliary equations
%
%  Called by functions: get_Gu.m, get_G_lambda.m, find_bifurcation_point.m 
%  arclength_newton_corrector.m, newton_corrector.m, switch_branch.m


if (p.sw.spcont==0)          % regular continuation 
    r=residual_pde(p,u);           % pde part
  if(p.nc.nq>0) 
      r=[r;p.fuha.qf(p,u)]; 
  end  % possibly add auxiliary functions
end

if p.sw.spcont==1 % bifurcation point continuation 
    ilam=p.nc.ilam; 
    u1=[u(1:p.nu);u(2*p.nu+1:2*p.nu+p.naux)]; % normal pde part
    r=residual_pde(p,u1); % pde part

    if(p.nc.nq>0) 
        rq=p.fuha.qf(p,u1); 
    else 
        rq=[]; 
    end
    
    G_lambda=get_G_lambda(p,u1,r);  Gu=get_Gu(p,u1,[r;rq]); 
    psi=u(p.nu+1:2*p.nu); mu=u(2*p.nu+ilam(end)); % the (zero) eigenvalue/vector 
    r=r+mu*p.mat.M(1:p.nu,1:p.nu)*psi; rlin=Gu'*psi; 
    rnorm=norm([psi; 0*u(2*p.nu+p.naux+1:end)])^2-1; rcon=psi'*G_lambda; 
    r=[r; rlin(1:p.nu); rq; rlin(p.nu+1:end); rnorm; rcon];   % sort into residual 
end



end