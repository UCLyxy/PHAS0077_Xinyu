function active_u=transform_u_activeu(p,varargin)
% transform_u_activeu: select only the active auxiliary variables and append to upde
%
% USAGE:
%  au=transform_u_activeu(p) - without primary parameter
%  au=transform_u_activeu(p,1) - with primary parameter for bordered system
%
% See also tansform_activeu_u

if ~(isempty(varargin))
    u=varargin{1};
else u=p.u; 
end

if (nargin>2 && varargin{2}==1) % bordered system with primary par
    active_u=zeros(p.nu+p.nc.nq+1,1);  % au=[upde,actaux,lam] 
    active_u(p.nu+p.nc.nq+1)=u(p.nu+p.nc.ilam(1)); % put primary parameter at end
else
    active_u=zeros(p.nu+p.nc.nq,1); % au1=[upde,actaux] 
end
active_u(1:p.nu)=u(1:p.nu); % upde entries

for j=1:p.nc.nq 
    active_u(p.nu+j)=u(p.nu+p.nc.ilam(j+1)); 
end % select active aux. var 

