function [p,iok]=initial_step(p)
% initial_step: return p-structure and ok flag for (attempt of) initial
% continuation step done by incrementing Parameter and Newton loop.
% This yields the first 2 points on branch and first tangent (secant) vector.
%
% USAGE:
%   [p,iok]=initial_step(p)
% 
% INPUT:
%   *p: object oriented pde object before first continuation
%   *iok: convergence error flag with value 0 or 1
% 
% Called by numerical_continuation.m

iok1=0;
iok2=0;
[u0,res,iter,~,~,p]=newton_corrector(p,p.u);
ineg=-1; muv=[]; 
p.sol.meth='natural'; % natural continuation
if(res<p.nc.tol) 
    p.u=u0; p.sol.res=res; p.sol.iter=iter; 
    if (p.sw.spcalc>0 || p.sw.bifcheck==2) % spectral calc 
        r=calculate_residual(p,p.u); Gu=get_Gu(p,p.u,r); [ineg,muv]=vspcalc(Gu,p); 
    end

    if(p.sw.errcheck>0)
        p.sol.err=errcheck(p);
    end

    if(~isfield(p.sol,'ptype') || p.sol.ptype~=-2) 
        p.sol.ptype=-1;
    end % initial point, unless from swibra

    brout=[default_branch_data(p); p.fuha.outfu(p,p.u)]; % userfu to append to bif-branches  
    brplot=brout(length(default_branch_data(p))+p.plot.bpcmp);
    p.branch=[p.branch brout]; % put on branch 
    figure(p.plot.brfig); hold on; % plot point 

    if(any(ineg<=0))
        plot(get_para_lambda(p),brplot,'*'); 
    else 
        plot(get_para_lambda(p),brplot,'+'); 
    end 

    p.sol.ineg=ineg; p.sol.muv=muv;

    if p.file.smod~=0 
        p.fuha.savefu(p);  
    end % save to file 
    dss=0; 
    [p,~]=p.fuha.ufu(p,brout,dss);     
    iok1=1;
    p.file.count=p.file.count+1; 
else 
    fprintf('   - no convergence in zeroth step, lam=%g, res=%g\n',get_para_lambda(p),res);
end
u0=p.u; lam0=get_para_lambda(p); p.u(p.nu+p.nc.ilam(1))=lam0+p.sol.ds/10; 

[u,res,iter,~,~,p]=newton_corrector(p,p.u); 
if(res<p.nc.tol) 
    
    p.u=u; p.sol.res=res; p.sol.iter=iter; 
    if(p.sw.spcalc>0 || p.sw.bifcheck==2); r=calculate_residual(p,u); Gu=get_Gu(p,u,r); 
        [ineg,muv]=vspcalc(Gu,p); 
    end
    
    if(p.sw.errcheck>0) 
        p.sol.err=errcheck(p);
    end 
    brout=[default_branch_data(p); p.fuha.outfu(p,p.u)]; % userfu to append to bif-branches  
    brplot=brout(length(default_branch_data(p))+p.plot.bpcmp);
    p.branch=[p.branch brout]; % put on branch 
    p.sol.ptype=0; % normal point
    figure(p.plot.brfig); hold on; % plot point 

    if(any(ineg<=0)) 
        plot(get_para_lambda(p),brplot,'*'); 
    else 
        plot(get_para_lambda(p),brplot,'+'); 
    end

    p.sol.ineg=ineg; p.sol.muv=muv; 
    tau=(transform_u_activeu(p,p.u,1)-transform_u_activeu(p,u0,1))/(get_para_lambda(p)-lam0);
    p.tau=tau/calculate_xi_norm(tau,p.sol.xi,p.nc.nq,p.sol.xiq); 
    dss=p.sol.ds/10; [p,~]=p.fuha.ufu(p,brout,dss);    
    iok2=1;p.file.count=p.file.count+1;
else 
    fprintf('   - no convergence in first step, lam=%g, res=%g\n',get_para_lambda(p),res);
end
iok=iok1*iok2; p.sol.restart=0; 


