function Gu=get_Gu(p,u,r)
% get_Gu: finite-difference computing,calculate partial derivatives of the finite-element approximation
% function G(u,lambda) (right hand side of equation) with respect to u
% (see dissertation equation(3.3))
%
%  USAGE:
%  Gu=get_Gu(p,u,r)
%
%  INPUTS:
%   *p: object oriented pde object
%   *u: solution vector of pde system(u=[u1;u2;par])
%   *r: residual computed by residual_pde.m
%
%  RETURNS:
%   *Gu: partial derivatives of the finite-element approximation function G(u,lambda) 
%        (right hand side of equation) with respect to u
%  
%  Reference: See dissertaion: total derivative equation (3.3)
%
% pde-part and auxiliart part separated.

Gu=p.fuha.sGjac(p,u);
if p.nc.nq>0 % derivatives w.r.t active auxiliary vars by finite diff.
    Gqw=zeros(p.nu+p.nc.nq,p.nc.nq); % (nu+nq) x nq
    qu=zeros(p.nc.nq,p.nu); 
    for k=2:p.nc.nq+1 % derivative wrt parameters, skip primary parameter
      r1=calculate_residual(p,u+p.nc.del*ej(p.nu+p.nc.ilam(k),length(u))); 
      Gqw(:,k-1)=(r1-r)/p.nc.del;
    end
    rq=r(p.nu+1:end); 
    if p.sw.qjac==1; qu=p.fuha.qfder(p,u); % analytical q_u 
    else % numerical q_u 
        for j=1:p.nu 
            r1=p.fuha.qf(p,u+p.nc.del*ej(j,length(u))); 
            qu(:,j)=(r1-rq)/p.nc.del; 
        end
    end
    Gu=[[Gu; qu] Gqw];
end


return;