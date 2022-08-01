function p=find_bifurcation_point(p,varargin)
% find_bifurcation_point: find and locate bifurcation points of p structure at current branch
%
% USAGE:
%  p=find_bifurcation_point(p)     - find one bifurcation using at most p.nc.nsteps
%  p=find_bifurcation_point(p,num_bif) - try to find num_bif bifurcation points using at most p.nc.nsteps
% 
% INPUTS:
%  *p: object oriented pde object
%  *num_bif: number of bifurcation points which user wants to find and
%  locate on the current branch of p 
%
% *Remark: this function should only be called on user command line
% interface to help user find and plot bifurcation points 

num_bif_point=1;

if nargin>1
    num_bif_point=varargin{1};
end

ds=p.sol.ds;

% if num_continuation has not yet been run make just one small step
if isfield(p.fuha, 'innerlss')
    lsss=p.fuha.innerlss; 
    p.fuha.innerlss=@lss; 
end 

if p.sol.restart==1
    p.file.pnamesw=0;
    p.sol.ds=1e-10; 
    try % this try-catch is when user does not want labels overwritten (set_file_name).
        p=numerical_continuation(p,1); p.file.count=p.file.count+1; 
    catch exception
        fprintf(exception.message); 
        return; 
    end
    p.sol.ds=ds;
end

% saving parameters, which will change locally to store back below
smod=p.file.smod;
bif_check=p.sw.bifcheck; 
fold_check=p.sw.foldcheck; 
pnamesw=p.file.pnamesw;
headfu=p.fuha.headfu; 
ufu=p.fuha.ufu; 
timesw=p.time.timesw; 
pmod=p.plot.pmod;
if p.sw.spcalc==0 % calculate ineg here since switched off in p 
   r=calculate_residual(p,p.u); Gu=get_Gu(p,p.u,r); 
   [p.sol.ineg,p.sol.muv]=spcalc(Gu,p); 
end
ineg0=p.sol.ineg; % number of negative eigenvalues of starting point
fprintf('Start find bifurcation points, %i negative EigenVals at lambda=%g\n', p.sol.ineg, get_para_lambda(p));
n=0; k=0; st=0; ft=p;
while (k<num_bif_point && st==0) % loop over nbp desired bif-points 
  k=k+1; af=0; p.file.smod=0; p.file.pnamesw=0; p.fuha.headfu=@findbifheadfu; % no saving/output
  p.fuha.ufu=@findbiffu; p.sw.bifcheck=0; p.sw.foldcheck=0; % the bifpoint will be calculated at the end 
  p.time.timesw=0; p.plot.pmod=0; p.sol.ds=ds; 
  while (st==0 && abs(p.sol.ds)>p.nc.dsminbis) % find next bifpoint
    if(abs(p.sol.ds)<p.nc.dsmin); st=1; fprintf('  ds<dsmin, stopping\n');end
    ft=p; ft=numerical_continuation(ft,1); n=n+1; 
    if(get_para_lambda(ft)<p.nc.lammin); st=2; fprintf('  lam<lammin, stopping\n');end
    if(get_para_lambda(ft)>p.nc.lammax); st=2; fprintf('  lam>lammax, stopping\n');end
    if(n>p.nc.nsteps); st=2; fprintf('  n>p.nc.nsteps, stopping\n');end
    if st<2
      if p.sw.spcalc==0 % need to calc ineg here since switched off in p 
      % of course somewhat inefficient, but we want to keep p.branch clean
       r=calculate_residual(ft,ft.u); Gu=get_Gu(ft,ft.u,r); [ft.sol.ineg,ft.sol.muv]=spcalc(Gu,ft); 
      end
      if(p.sw.verb>0); fprintf('   lam=%s, ineg=%i\n',printcon(get_para_lambda(ft)),ft.sol.ineg); end    
      if ft.sol.ineg==ineg0 % index didn't change, next step! 
        p=ft;
        if af==1; p.sol.ds=p.sol.ds/2; end % index-change already found, decrease ds 
        af=0; % assume that next solution will have no large index-change
      else
        if abs(ft.sol.ineg-ineg0)>1 % index-change by more than 1 found
           af=1; if abs(p.sol.ds)>p.nc.dsminbis; p.sol.ds=p.sol.ds/2; end; % decrease ds 
        else 
           st=1; % point found! 
        end
      end
    end
  end   % while find next bifpoint
  if (p.sol.deta==0) p.sol.deta=get_determainant(p); end % p.sol.deta not yet calculated  
  % now bifurcation (unless index change due to fold) should be found 
  % with current p.sol.ds, first restore parameters 
  ineg0=ft.sol.ineg; 
  p.file.smod=smod; p.sw.bifcheck=bif_check; p.sw.foldcheck=fold_check; 
  p.fuha.headfu=headfu; p.fuha.ufu=ufu; 
  p.time.timesw=timesw; p.plot.pmod=pmod; 
  if st==1; p.file.count=p.file.count+1; p=numerical_continuation(p,1); st=0; end % find bif.point using cont
end % loop over nb (desired) bifpoints
% set remaining parameters back to original
p.sw.bifcheck=bif_check; p.file.pnamesw=pnamesw; p.sol.ds=ds; 
if isfield(p.fuha, 'innerlss'); p.fuha.innerlss=lsss; end  
end 

% inner functions
function [p,cstop]=findbiffu(p,~,~)
    p=para_lambda_test(p); % test if a desired lambda value has been passed
    cstop=0;
end

function findbifheadfu(~)
end
