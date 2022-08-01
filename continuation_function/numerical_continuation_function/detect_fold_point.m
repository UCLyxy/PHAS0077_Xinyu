function p=detect_fold_point(p,u1,tau1)
% detect_fold_point: called from num_continuation to check for 
% the change of sign(det(Gu)) as well as the Jacobian of G.
%   *If fold point is detected, locate by bisection method and save it.
%
% USAGE:
%       p=detect_fold_point(p,u1,tau1)
%
% 
% p.u, p.sol.deta, ineg0   at previous continuation point   
% u1, tau1, Gu, Glam       at new point.
%
% Called by function: numerical_continuation.m 


ds=p.sol.ds; % continuation step size

if(p.sol.ptype<0 || p.sol.restart>0 || abs(sign(tau1(p.nu+p.nc.nq+1))-sign(p.tau(p.nu+p.nc.nq+1)))<2) 
    % restart, or no fold, do nothing 
    return; 
end

% Fold detected
if(p.sw.verb>0)
    fprintf('   detect fold point between %g and %g\n',get_para_lambda(p),get_para_lambda(p,u1)); 
end 

% Use bisection method to locate
u0=p.u;
tau0=p.tau; 
ds=ds/2; 
bisecc=0; 
% start bisection  
while (bisecc<p.nc.bisecmax && abs(ds)>p.nc.dsminbis)  
   au0=transform_u_activeu(p,u0,1); au1=transform_u_activeu(p,u1,1); % get active variables
   if p.sw.bifloc==0 % tangent predictor 
      aun=au0+ds*tau0; %select entries of tau0 into full u = [upde;par]
   else if p.sw.bifloc==1 % secant predictor 
          aun=au0+0.5*(au1-au0); 
        else % quadratic predictor
          aun=0.25*(3*au0+au1)+0.5*ds*tau0; 
       end
   end
   un=transform_activeu_u(p,aun,1); u0=transform_activeu_u(p,au0,1); % merge with auxiliary variables

   if(p.sw.verb>1)
       fprintf('   checking lam=%g ...',aun(p.nu+p.nc.nq+1));
   end
   
   if(p.sw.para==0 || (p.sw.para==1 && abs(tau0(p.nu+p.nc.nq+1))>p.nc.lamdtol)) % fixed lam corrector
      [un,res,~,Gu,Glam]=newton_corrector(p,un); 
   else
      [un,res,~,Gu,Glam]=arclength_newton_corrector(p,un,ds,u0); % arclength corrector
   end

   if(res<2*p.nc.tol) % If step was (rather) OK, then 
     arclength_jac_matrix=gen_ext_mat(p,Gu,Glam,tau0,p.sol.xi,p.sol.xiq); 
     taun=arclength_jac_matrix\[zeros(p.nu+p.nc.nq,1);1]; 
     taun=taun/calculate_xi_norm(taun,p.sol.xi,p.nc.nq,p.sol.xiq); % tau at new point

     if(p.sw.verb>1) 
         fprintf('   ok, sign(lamd)=%i\n',sign(taun(p.nu+p.nc.nq+1))); 
     end 

     if sign(taun(p.nu+p.nc.nq+1))==sign(tau0(p.nu+p.nc.nq+1)) % u_n still before bif, reset u0, keep u1; 
         ds=ds/2; 
     else  % u_n is behind bif., reset u1 to u0, i.e., flip direction:
        u1=u0; 
        ds=-ds/2; 
     end
     u0=un; tau0=taun; 
   else
     if(p.sw.verb>1)
         fprintf('No convergence, localization might be poor ...\n');
     end
     ds=3*ds/4; % might have hit sing point, try a smaller step in predictor
   end % if res<2*p.nc.tol 
   bisecc=bisecc+1;
end % localization done 
%------------------------------------------------------- postprocessing!
us=p.u; taus=p.tau; tys=p.sol.ptype; % remember values behind fold 
p.u=u0; p.tau=tau0; % and store current values for save 
p.sol.ptype=2; % fold point

if(p.sw.errcheck>0)
    p.sol.err=errcheck(p);
end 

brout=[default_branch_data(p); p.fuha.outfu(p,u0)]; % userfu to append to bif-branches  
brplot=brout(length(default_branch_data(p))+p.plot.bpcmp);
p.branch=[p.branch brout]; p.fuha.savefu(p); % put on branch and save   
fname=[p.file.fpname,sprintf('%i',p.file.fcount),'.mat']; %save(fname,'p');    
fprintf('%4i %s (FP, saved to %s) bisection steps %i, last ds %g\n',p.file.count,printcon(get_para_lambda(p)),...
    fname,bisecc,ds);
figure(p.plot.brfig); plot(get_para_lambda(p),brplot,'x');
if((p.file.count>0) && (mod(p.file.count,p.file.smod)==0)) % save to file 
   fname=[p.file.pname,sprintf('%i',p.file.count),'.mat']; % save(fname,'p'); 
   p.fuha.savefu(p,0); fprintf('   also saved to %s\n',fname);
end 
p.u=us; p.tau=taus; p.sol.ptype=tys; % restore values behind bif 
p.file.fcount=p.file.fcount+1; p.file.count=p.file.count+1; 


