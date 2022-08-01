function [p,cstop]=stanufu(p,brout,ds)
% STANUFU: standard "user function" called after each continuation step
% USAGE:
%  [p,cstop]=stanufu(p,brout,ds)
%
% INPUTS:
%   *brout:branch data from internal calculations
%   *ds:current stepsize
% Returns cstop flag for error post-processing
%
% Called by numerical_continuation.m

p=para_lambda_test(p); % test if a desired lambda value has been passed
brplot=brout(length(default_branch_data(p))+p.plot.bpcmp); %y-axis value in bif-figure
fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
    p.file.count,printcon(get_para_lambda(p)),printcon(brplot),p.sol.res,p.sol.iter, ...
    p.sol.meth,printcon(ds));

if(p.sw.errcheck>0) 
    fprintf('%5.2e ',p.sol.err);
end

if(p.sw.spcalc==1)
    fprintf(' %2i ', p.sol.ineg); 
end

npr=length(p.sw.bprint); 

for i=1:npr
    fprintf('%s ',printcon(brout(length(default_branch_data(p))+p.sw.bprint(i)))); 
end
% put anything else here
fprintf('\n');
cstop=0;
if(get_para_lambda(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',get_para_lambda(p),p.nc.lammin); cstop=1;
end
if(get_para_lambda(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',get_para_lambda(p),p.nc.lammax); cstop=1;
end
end
