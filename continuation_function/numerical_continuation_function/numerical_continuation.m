function p=numerical_continuation(p,num_continuation)
% numerical_continuation: main continuation routine 
%
% USAGE:
%   p=numerical_continuation(p)   :  do p.nc.nsteps continuation steps
%   p=numerical_continuation(p,n) :  do n continuation steps
%
% INPUT:
%   *p: object oriented pde structure
%   *num_continuation: the number of continuation steps that user would
%   like to take
% 
% MAIN ALGORITHM:
%   *compute predictor: (u1,lam1)=(u0,lam0)+ds*tau0
%   *apply Newton corrector:use newton_corrector.m for natural continuation
%    and arclength_newton_corrector.m for aclength continuation(fold point)
%   *compute new tangent vector tau1
%
% For details, see dissertation Section 3.3
%
% Remark: Initially, make continuation with stepsize p.sol.ds, starting from p.u, direction tau
% save the initial p object in case of failure

back_up_p=p; 

%Check whether the mass matrix of fe is defined or not
if p.sw.sfem~=0 
    if ((p.sw.spcalc~=0 && isempty(p.mat.M)))  
    fprintf('\nMass matrix is not defined, call p=setfemops(p)\n');
    p=setfemops(p); 
    end
end

% set variable names to the file
if(p.file.pnamesw==1) 
    [p,bool]=set_file_name(p,inputname(1)); 
    if bool~=1; return; end
end

% set labels to axes
x_length=0;

if(isfield(p.plot,'auxdict')) 
    x_length=length(p.plot.auxdict); 
end

if(p.nc.ilam(1)<=x_length) 
    x_axis=p.plot.auxdict{p.nc.ilam(1)};  
else 
    x_axis=['parameter=' mat2str(p.nc.ilam(1))];
end

figure(p.plot.brfig);
xlabel(char(x_axis)); 
y_axis=mat2str(p.plot.bpcmp);

if(p.plot.bpcmp>0)
  if(p.plot.bpcmp<=x_length) 
      y_axis=p.plot.auxdict{p.plot.bpcmp};  
  else 
      y_axis=['user branch comp. ' y_axis]; 
  end
else
    y_axis='L2-norm';
end
ylabel(char(y_axis));

%init timing and record
p.time.newton=0; 
p.time.totst=0; 
p.time.spec=0; 
totime=tic; 
max_steps=p.nc.nsteps;

% define continuation step
if nargin>1
    max_steps=num_continuation; 
end

% count iteration
num_iteration=0;
% record last negative eigenvalues
ineg0=p.sol.ineg; 

% output header created by user
p.fuha.headfu(p);

% loop for numerical continuation
while num_iteration<max_steps
    % initial iteration
   if (length(p.tau)~=p.nu+p.nc.nq+1)||(p.sol.restart>0)

      if(p.file.count>0) 
          p.file.count=p.file.count-1;
      end

      [p,iok]=initial_step(p);

      if iok==0
          p=back_up_p; 
          return; 
      end

      ineg0=p.sol.ineg; 
      plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); 
      num_iteration=num_iteration+1;

      if num_iteration>=max_steps
          p.file.count=p.file.count-1; 
          p.fuha.savefu(p); 
          return; 
      end
   end
    
   step_time=tic;
   num_iter=0;
   residual=10*p.nc.tol;
   % step_bool=1 after successful step 
   step_bool=0;
   % save ds of current iteration 
   ds_vec=p.sol.ds; 
   while step_bool==0     % iteration to find next solution point
      para_lambda=p.tau(p.nu+p.nc.nq+1);
      % define the predictor
      active_u1=transform_u_activeu(p,p.u,1)+p.sol.ds*p.tau;
      u1=transform_activeu_u(p,active_u1,1); 
      newton_time=tic;

      if(p.sw.para==0 || (p.sw.para==1 && abs(para_lambda)>p.nc.lamdtol)) 
          % natural lambda corrector
          [u1,residual,num_iter,Gu,Glam,p]=newton_corrector(p,u1); p.sol.meth='natural';

      else 
          % arclength-corrector  
          [u1,residual,num_iter,Gu,Glam,p]=arclength_newton_corrector(p,u1,p.sol.ds); 
          p.sol.meth='arclength';  
      end
      % calculate newton-iteration time
      p.time.newton=p.time.newton+toc(newton_time);
      % record step size
      ds_vec=p.sol.ds;       
      [p,step_bool,u1,residual,num_iter,Gu,Glam]=stepsize_control_continuation(p,u1,residual,num_iter,Gu,Glam,ds_vec); 

   end              
   p.sol.ptype=0; % so far normal point
   if p.sw.spcalc>0 % calculate eigenvalues 
       sptime=tic; ineg0=p.sol.ineg; [p.sol.ineg,p.sol.muv]=vspcalc(Gu,p); 
       p.time.spec=p.time.spec+toc(sptime); % spectral-time, accumulated 
   end
   secpred=0; try secpred=p.sw.secpred; catch; end 
   if secpred==1 % secant instead of tangent 
     u_a2=transform_u_activeu(p,u1,1);
     u_a1=transform_u_activeu(p,p.u,1); 
     tau_pre=sign(p.sol.ds)*(u_a2-u_a1); 
   else  % form extended matrix and compute new tangent
     amat=gen_ext_mat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); 
     [tau_pre,p]=p.fuha.blss(amat,[zeros(p.nu+p.nc.nq,1);1],p); 
   end
   tau_pre=tau_pre/calculate_xi_norm(tau_pre,p.sol.xi,p.nc.nq,p.sol.xiq);
   file_count=p.file.count;    
   if ~isfield(p.sw,'abs'); p.sw.abs=0; end % abs=1 if directly after a bifurcation 
   if p.sw.bifcheck>0 && p.sw.abs~=1  % check for bifurcation, unless 1st step after bif
     next_ds=p.sol.ds; p.sol.ds=ds_vec; 
     biftime=tic; 
     [p,bif_bool]=detect_bifurcation_point(p,u1,tau_pre,Gu,Glam,ineg0); 
     ineg0=p.sol.ineg; 
     p.time.bif=p.time.bif+toc(biftime); 
     if bif_bool; if p.sw.cdbb==1; u1=p.u; p.sol.ds=next_ds/2; end
     else p.sol.ds=next_ds; 
     end
   end

   if(p.sw.foldcheck>0 && abs(sign(tau_pre(p.nu+p.nc.nq+1))-sign(para_lambda))>1 &&p.sw.abs~=1)  
     next_ds=p.sol.ds; p.sol.ds=ds_vec; 
     p=detect_fold_point(p,u1,tau_pre);
     p.sol.ds=next_ds; 
   end   
    
   num_iteration=num_iteration+(p.file.count-file_count);
   p.sol.restart=0; p.sw.abs=0; p.u=u1; p.tau=tau_pre;     
   p.sol.res=residual; p.sol.iter=num_iter; p.sol.lamd=para_lambda; 
   if(p.sw.errcheck>0); p.sol.err=errcheck(p);  
     if(p.sol.err>p.nc.errbound && p.nc.errbound>0)                            
       if(p.sw.errcheck==1  || p.sw.bcper~=0) 
         fprintf('   - err.est.=%g>errbound=%g. Consider mesh-refinement.\n', p.sol.err, p.nc.errbound);          
       end 
       if(p.sw.errcheck==2 && p.sw.bcper==0)
           fprintf('   - err.est.=%g>errbound=%g. Adapting mesh\n', p.sol.err, p.nc.errbound); 
           p=meshadac(p,'eb',p.nc.errbound); 
       end 
       if(p.sw.errcheck>2) && p.sw.bcper==0; p=meshref(p,'eb',p.nc.errbound); end % refine mesh 
     end 
   end   
   
   if (mod(p.file.count,p.nc.amod)==0 && p.file.count~=0)
       fprintf('   - adap mesh\n'); 
     p=meshadac(p,'eb',0);  
   end    
   brout=[default_branch_data(p); p.fuha.outfu(p,p.u)];          % userfu to append to bif-branches  
   brplot=brout(length(default_branch_data(p))+p.plot.bpcmp);    %y-axis value in bif-figure
   p.branch=[p.branch brout];                       % put on branch 
   figure(p.plot.brfig); hold on;                   % plot point 
   if p.sol.ineg<=0; plot(get_para_lambda(p),real(brplot),'*b'); drawnow; 
   else plot(get_para_lambda(p),real(brplot),'+b'); drawnow; end 
   p.time.totst=p.time.totst+toc(step_time);            % total step time  
   [p,cstop]=p.fuha.ufu(p,brout,ds_vec);                
   if(p.file.count>0 && mod(p.file.count,p.file.smod)==0) % save in file 
      p.fuha.savefu(p); end
   if(mod(p.file.count,p.plot.pmod)==0) 
       plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle);end % plot solution       
   p.file.count=p.file.count+1; num_iteration=num_iteration+1; 
   if p.file.count>p.nc.ntot; cstop=1; end 
   if isfield(p.mat,'prec')
       try
   if(mod(p.file.count-1,p.ilup.ilun)==0) % force update of prec
     disp('new prec');
     p.mat.prec=AMGdelete(p.mat.prec); p.mat=rmfield(p.mat,'prec');   
     [p.mat.prec,p.amgopt]=AMGfactor(Gu,p.amgopt);
   end
       end
   end
   if cstop==1; break; end                           
end % while is<msteps
p.file.count=p.file.count-1; % some postprocessing, i.e., save the last point 
if(mod(p.file.count,p.file.smod)~=0 && p.file.smod~=0); p.fuha.savefu(p); end % save last point with adjusted counter 
p.file.count=p.file.count+1; p.time.tot=toc(totime); 
if(p.time.timesw>0); fprintf('Timing: total=%g, av.step=%g, av.Newton=%g, av.spcalc=%g\n',...
      p.time.tot,p.time.totst/num_iteration,p.time.newton/num_iteration,p.time.spec/num_iteration);
end
