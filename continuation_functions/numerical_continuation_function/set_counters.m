function p=set_counters(p)
% set_counters: reset counters in p and clear p.branch
% USAGE:
%  p=set_counters(p)
%
% Called by function: switch_branch.m

p.branch=[]; 
p.file.count=0; 
p.file.bcount=1; 
p.file.fcount=1; 

