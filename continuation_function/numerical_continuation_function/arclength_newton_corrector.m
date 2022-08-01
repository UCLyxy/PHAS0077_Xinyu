function [u,res,num_iter,Gu,G_lambda,p]=arclength_newton_corrector(p,u1,ds,varargin)
% arclength_newton_corrector: newton method applied for the case: arclength 
% continuation(fold point) corrector
%
% USAGE:
%   [u,res,num_iter,Gu,G_lambda]=arclength_newton_corrector(p,u1,ds)   : use p.u and u1 for initial arclength step
%   [u,res,num_iter,Gu,G_lambda]=arclength_newton_corrector(p,u1,ds,u) : use u and u1 for initial arclength step
%
% INPUTS:
% * u1 is initial guess, ds initial stepsize usually taken from p-structure.
%
% OUTPUTS:
% * Returns u, res=residual, num_iter=number of iterations, Gu,G_lambda=derivatives.
%
% * Remark: returns full vector u (including auxiliary variables such as u2)! 
% * Critical parameter settings:  p.nc.tol, p.nc.imax, p.sw.newt, p.fuha.blss
% (See stanpara.m for detail information of parameter setting)
%
% Reference: dissertation(Section 3.3 corrector)
% Called by functions: numerical_continuation.m

if ~isempty(varargin)
    active_u0=transform_u_activeu(p,varargin{1},1);
else 
    active_u0=transform_u_activeu(p,p.u,1); % last point on branch
end

num_iter=0;
u=u1;
alpha=1;

try almin=p.nc.almine; 
catch 
    almin=0.2; 
end % minimal damping, try-catch for backward-comp.

r=calculate_residual(p,u);
res0=norm(r,p.sw.norm); 
res=res0; % the residual 

[Gu,G_lambda]=get_info_derivative(p,u,r); % the derivatives 

if(res<p.nc.tol)
    return; 
end % initial res. small, do nothing!

%---------------------------------------------------------- start newton method 
step_bool=1; % step_bool=1 indicates that step of size alpha decreased residual!
             % step_bool=0 indicates that step of size alpha increased
             % residual

while(abs(res0)>p.nc.tol && num_iter<p.nc.imax && step_bool) % stopping criteria

  arclength_jac_matrix=gen_ext_mat(p,Gu,G_lambda,p.tau,p.sol.xi,p.sol.xiq);
  % see equation (3.13) of the dissertation
  active_u=transform_u_activeu(p,u1,1); 
  aud=active_u-active_u0; 

  if(p.nc.nq>0)
    p1=p.tau'*[p.sol.xi*aud(1:p.nu);p.sol.xiq*aud(p.nu+1:p.nu+p.nc.nq);(1-(p.sol.xi+p.sol.xiq)/2)*aud(p.nu+p.nc.nq+1)]-ds;
  else
    p1=p.tau'*[p.sol.xi*aud(1:p.nu);(1-p.sol.xi)*aud(p.nu+p.nc.nq+1)]-ds;
  end  

  [upd,p]=p.fuha.blss(arclength_jac_matrix,[r;p1],p);
  step_bool=0; 

  while(step_bool==0 && alpha>almin) 
    au1=active_u-alpha*upd;
    u1=transform_activeu_u(p,au1,1); 
    r=calculate_residual(p,u1); 
    res=norm(r,p.sw.norm); 
    if(res<res0) % good step 
       if(res<res0/2 && alpha<1) 
           alpha=alpha*2; 
       end % very good step, possibly increase alpha 

       step_bool=1; 
       u=u1; 
       res0=res;

       if(p.sw.newt==0) 
           [Gu,G_lambda]=get_info_derivative(p,u,r); 
       end % full Newton, get new derivatives 
    else 
        alpha=alpha/2; % bad step, try smaller alpha
    end % res<res0
  end % while stepok==0
  num_iter=num_iter+1;
end % while res<p.nc.tol && stepok 
% ----------------------------------------------------------postprocessing
if(p.sw.newt>0)  % if chord, then now get derivatives at new point! 
    [Gu,G_lambda]=get_info_derivative(p,u,r);
end

if(p.sw.verb>1)
    if(alpha<1) 
        fprintf('\narclength_newton_corrector: damp alpha=%g, res=%g, ds=%g\n',...
            alpha,res,ds); 
    end 
end % inform user if damping was used 
res=res0; 
