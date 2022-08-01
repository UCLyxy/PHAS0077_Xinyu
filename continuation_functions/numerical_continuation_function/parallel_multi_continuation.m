function p=parallel_multi_continuation(p,varargin)
% parallel_multi_continuation: parallel multi continuation 
%
%  p=parallel_multi_continuation(p)
% do p.nstep multisteps of (initial) size i*p.ds, i=1:p.pm.mst,
% Supplementary files for numerical_continuation

if nargin>1; p.nc.ntot=varargin{1}; end 
cstop=0; is=0; p0=p; 
m=p.pm.mst; scount=p.file.count; % save counter for timing

if(length(p.nc.ilam)~=p.nc.nq+1)
    fprintf('\nLength of p.nc.ilam does not match equations plus primary parameter!\n');
    fprintf('\nParameter setting error!\n');
return; 
end

if( (p.sw.spcalc~=0 && isempty(p.mat.M)) || (p.sw.sfem~=0 && (isempty(p.mat.M)) )) 
    fprintf('\nNo defined mass matrix  -- calling function p=setfemops(p)\n');
    p=setfemops(p);
end

if(p.file.pnamesw==1) % check if dir exists
    dum1=p.file.count; dum2=p.branch; dum3=p.file.bcount; 
    [p,ok]=set_file_name(p,inputname(1)); 
    p.file.count=dum1; p.branch=dum2; p.file.bcount=dum3; 
    if ok~=1 
        return;
    end
end

figure(p.plot.brfig); 
p.time.newton=0; p.time.st=0; % timers for parallel_multi_newton and multisteps! 
totime=tic; % total time timer 
%gcp; %opens cores for parallel computing, if not open 
p.fuha.headfu(p); 
for k=1:p.nc.nsteps % the outer steps 
  stime=tic; 

  if ((length(p.tau)~=p.nu+p.nc.nq+1)||(p.sol.restart>0)) % initial step 
      [p,iok]=inistep(p);
      if(iok==0)
          return;
      end 
      plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle);
      title('');
      is=is+2;
  end

  ineg0=p.sol.ineg; % use last ineg for bifdetec consistency check       

  p.pm.imax=p0.pm.imax;   % set p.pm.imax to orignal one 
  foso=0; % foso=found solution   
  while foso==0  % change, e.g., ds, if foso=0 in parallel_multi_newton 
     ntime=tic;  % timer for mst newton-loops! 
     
     [uold,tauold,inegp,resp,dsp,pp]=parallel_multi_newton(p,p0); % main call *********************
     
     p.time.newton=p.time.newton+toc(ntime);
     % time for p.pm.mst newtonloops 
     for i=1:p.pm.mst
         if (resp(i)<p.nc.tol)
             foso=foso+1; 
         end
     end  % get # of found solutions 

     if (foso==0 && abs(p.sol.ds)<p.nc.dsmin*(p.pm.mst+1) && p.nc.imax>p.pm.imax) 
         p.pm.imax=p.pm.imax+1; 
     end

     if (foso==0 && abs(p.sol.ds)>=p.nc.dsmin*(p.pm.mst+1))
         p.sol.ds=p.sol.ds/(p.pm.mst+1); 
     end
 
     if (foso==0 && abs(p.sol.ds)<p.nc.dsmin*(p.pm.mst+1) && p.nc.imax==p.pm.imax) 
         disp('no solution found'); 
         foso=-1; 
     end
   end % while foso 

   p.nc.imax=p0.nc.imax; % "postprocessing", check bif, plot, calc. tangent
   p.sol.ptype=0; % so far normal point 
   msucc=0; lamd=p.tau(p.nu+p.nc.nq+1); 
   %ineg0, inegp, pause
   for i=1:m % check bif, plot, calc. tangent 
    if resp(i)<p.nc.tol        
       msucc=msucc+1; p.sol.restart=0; % step accepted!     
       Gu=pp(i).Gu; Glam=pp(i).Glam; p.sol.iter=pp(i).iter; ds=dsp(i); 
       p.sol.res=resp(i); lam1=get_para_lambda(pp(i)); u1=pp(i).u; 
       p.sol.meth=pp(i).meth; p.sol.ineg=inegp(i); % HU, 17-09
       p.sol.err=pp(i).sol.err; % HU 
       arclength_jac_matrix=gen_ext_mat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); 
       tau=p.fuha.blss(arclength_jac_matrix,[zeros(p.nu+p.nc.nq,1);1],p); 
      
        tau=tau/calculate_xi_norm(tau,p.sol.xi,p.nc.nq,p.sol.xiq);

        if ~isfield(p.sw,'abs')
            p.sw.abs=0; 
        end

       if(p.sw.bifcheck>0 && p.sw.abs~=1) % check for bif 
         p.u=uold;
         p.tau=tauold; 
         if(i>1); ineg0=inegp(i-1,:);
         end 
         p=detect_bifurcation_point(p,u1,tau,Gu,Glam,ineg0); 
       end

       if(p.sw.foldcheck>0 && abs(sign(tau(p.nu+p.nc.nq+1))-sign(lamd))>1)  % fold via sgn(lamd)
        newds=p.sol.ds; p.sol.ds=ds; % use ds from BEFORE last stepsize_control_continuation 
        p=detect_fold_point(p,u1,tau); p.sol.ds=newds; 
       end
       p.u=u1; p.tau=tau; p.sol.lamd=p.tau(p.nu+1); lamd=p.sol.lamd;     
       brout=[default_branch_data(p); p.fuha.outfu(p,p.u)]; % userfu to append to bif-branches  
       brplot=brout(length(default_branch_data(p))+p.plot.bpcmp);
       p.branch=[p.branch brout]; % put on branch 
       figure(p.plot.brfig); hold on; % plot point

       if(p.sol.ineg<=0)
           plot(lam1,brplot,'*');drawnow;
       else 
           plot(lam1,brplot,'+');drawnow; 
       end

       if(p.file.count>0 && mod(p.file.count,p.file.smod)==0) % save to file 
             p.fuha.savefu(p);
       end

       if(mod(p.file.count,p.plot.pmod)==0) 
           plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); 
           title(''); 
           axis tight;  
       end % plot sol

       uold=u1; 
       tauold=tau; %for bifdetec

       [p,cstop]=p.fuha.ufu(p,brout,ds);

       p.file.count=p.file.count+1; 
       is=is+1; 
       if p.file.count>p.nc.ntot
           cstop=1; 
       end

       if cstop==1; break; end
    end %res<p.tol 
    p.nc.imax=p0.nc.imax; %because is used in Newton loop and bifdetect
  end % i=1:m in postprocessing 

  if cstop==1 
      break; 
  end

  if (foso==p.pm.mst && abs(p.sol.ds)<=p.nc.dsmax/p.nc.dsincfac) 
     p.sol.ds=p.nc.dsincfac*p.sol.ds; 
 
  end

  if foso==-1 
      break; 
  end  


  cstime=toc(stime); p.time.totst=p.time.totst+cstime;

  if p.time.timesw>0 
      fprintf('time for %2i successful steps: %g\n',msucc,cstime);
  end
  
end % nstep loop
p.file.count=p.file.count-1; 
if(mod(p.file.count,p.file.smod)~=0&&p.file.smod~=0); p.fuha.savefu(p); end % save last point with adjusted counter
p.file.count=p.file.count+1; 
p.time.tot=toc(totime); 
if p.time.timesw>0 
    fprintf('Total time=%g, av.time per succ.step=%g, av.newton-time=%g\n', ... 
        p.time.tot, p.time.tot/(p.file.count-scount),p.time.newton/(p.file.count-scount)); 
end
