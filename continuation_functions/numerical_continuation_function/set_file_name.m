function [p,ok]=set_file_name(p,varargin)
% set_file_name: set problem directory name
%
% USAGE:
%  [p,ok]=set_file_name(p)      - use variable name as directory name
%  [p,ok]=set_file_name(p,dir)  - use "dir"
%
% *ok:flag for error catching
%
% Remark: user could call this function in command line interface to give a
% new name of problem directory
%
% Called by switch_branch.m, numerical_continuation.m

ok=1;

if isempty(varargin)
    dir=sprintf('%s',inputname(1)); 
else 
    dir=varargin{1}; 
end 

fprintf('Problem directory name: %s\n',dir);

if ~exist(dir,'dir') 
    mkdir(dir); 
    fprintf('creating directory %s\n',dir);
end

p.file.dir=dir; 
p.file.pname=[dir '/pt']; % solution point
p.file.bpname=[dir '/bpt']; % bifurcation point
p.file.fpname=[dir '/fpt']; % fold point