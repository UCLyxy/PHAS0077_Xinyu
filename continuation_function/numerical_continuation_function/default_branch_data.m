function branch_data=default_branch_data(p)
% default_branch_data: set default non-user-defined branch data 
%
%  USAGE:
%   branch_data=default_branch_data(p)
%
%  INPUTS:
%   *p: object oriented pde object
%
%  RETURNS:
%   *branch_data: default non-user-defined branch data
%   *branch_data=[p.file.count; p.sol.ptype; p.sol.ineg; get_para_lambda(p); 
%      p.sol.err; L2-norm 1st component]
%
%  Called by functions: detect_bifurcation_point.m, detect_fold_point.m
%                       initial_step.m, numerical_continuation.m

ineg=p.sol.ineg; % number of negative eigenvalues
ineg=max(ineg);
u=p.u; 
M=getM(p); 
n1=floor(p.nu/p.nc.neq); 
l2=sqrt(u(1:n1)'*M(1:n1,1:n1)*u(1:n1)); 
branch_data=[p.file.count; p.sol.ptype; ineg'; get_para_lambda(p); p.sol.err; l2]; 
p.bra=branch_data;  