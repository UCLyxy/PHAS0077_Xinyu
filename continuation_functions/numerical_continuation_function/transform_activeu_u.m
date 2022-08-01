function u=transform_activeu_u(p,active_u,varargin)
% transform_activeu_u: transform active u to u with active and inactive parts
%
% USAGE:
%   u=transform_activeu_u(p,au,varargin)
% upde and active auxiliary variables sorted into u
% 
% Cases:
% varargin=0 or none: without primary parameter (e.g. for lss)
% varargin=1: with primary parameter (e.g. for blss)
%
% See also transform_u_activeu


u=p.u;
u(1:p.nu)=active_u(1:p.nu);
% init and upde entries

for j=1:p.nc.nq
    u(p.nu+p.nc.ilam(j+1))=active_u(p.nu+j);
end

% aux except primary param
if (~isempty(varargin) && varargin{1}==1) 
    u(p.nu+p.nc.ilam(1))=active_u(p.nu+p.nc.nq+1); 
end
