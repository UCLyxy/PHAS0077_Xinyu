function out=stanbra(p,u)
% STANBRA: standard function for output beyond default_branch_data(p) to p.branch.
% here output takes the form: [pars, max(abs(u1)), min(abs(u1))
% higly problem dependent, take this as a template
%
upde=u(1:p.nu); 
np=p.nu/p.nc.neq; 
out=[u(p.nu+1:end); max(abs(upde(1:np))); min(abs(upde(1:np)))]; 