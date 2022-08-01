function [u,res,iter,Gu,G_lambda,p]=newton_corrector(p,u1)
% newton_corrector: newton method applied for the case: natural continuation corrector 
%  
%  USAGE:
%  [u,res,iter,Gu,Glam,p]=newton_corrector(p,u1)
%
% INPUTS:
%   *p: object oriented pde object
%   *u1: initial guess of pde, original solution part
% OUTPUTS
%   *Returns resulting u
%   *res: residual
%   *iter: number of iterations
%   *Gu,Glam: derivative information
%
% * Note: returns full vector u (including auxiliary variables)! 
% * Critical parameter settings:  p.nc.tol, p.nc.imax, p.sw.newt, p.fuha.blss
% (See stanpara.m for detail information of parameter setting)
%
%
% Reference: dissertation(Section 3.2)
% Called by functions: numerical_continuation.m


try almin=p.nc.almin;
catch 
    almin=0.2;
end % minimal damping, try-catch for backward-comp.

u=u1; 
alpha=1;
iter=0; 
r=calculate_residual(p,u); res0=norm(r,p.sw.norm); % (starting) residual 
Gu=get_Gu(p,u,r); % derivatives (needed for next tangent, thus always computed) 

if res0<p.nc.tol % starting residual already small, compute Glam and return 
   G_lambda=get_G_lambda(p,u,r); res=res0; return; 
end 
% now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual! 
while(res0>p.nc.tol && iter<p.nc.imax && stepok)  
  [upd,p]=p.fuha.lss(Gu,r,p); stepok=0; % iter, pause
  while(stepok==0 && alpha>almin)       
    au1=transform_u_activeu(p,u,1) ... % au1=[upde,actaux,lam] 
        -alpha*[upd;0]; % the newton step, no change in primary param. 
    u1=transform_activeu_u(p,au1,1); 
    r=calculate_residual(p,u1); res=norm(r,p.sw.norm);  % residual     
    if res<res0 % good step 
       if(res<res0/4 && alpha<1) % very good step, possibly increase alpha 
           alpha=2*alpha; 
       end 
       stepok=1; u=u1; res0=res; 
       if(p.sw.newt==0) % full Newton, get new derivatives 
          Gu=get_Gu(p,u,r);  % next Jacobian 
       end
    else 
        alpha=alpha/2; % bad step, try smaller al
    end % if res<res0
  end % while stepok==0 
  iter=iter+1; %alpha=1; 
end   % while res0>p.nc.tol && stepok 
% some postprocessing

if(p.sw.newt==0)
    G_lambda=get_G_lambda(p,u,r); % Newton, only update Glam 
else 
    [Gu,G_lambda]=get_info_derivative(p,u,r);
end % chord, update Gu,Glam 

if(p.sw.verb>1); if(alpha<1) % inform user about damping ...
        fprintf('newton_corrector: damp alpha=%g, res=%g\n', alpha,res0); end; end; 
res=res0; % return res of u0=best u
