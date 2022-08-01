function p=setfemops(p)
% SETFEMOPS: generate and store FEM operators 
% Uses full domain assembling, and transforms in case of periodic bc.
%
%  p=setfemops(p)
%

if p.sw.sfem<0
   if isfield(p.fuha,'setops'); p=p.fuha.setops(p); % possible function handle 
   else p=oosetfemops(p); % for backward compatibility 
   end
   return; 
end
try upde=p.mat.fill*p.u(1:p.nu); % set to full domain vector
catch; upde=zeros(p.nu,1); end 
% fold/branch point continuation: set neq to original value:
if (p.sw.spcont~=0) 
    neq=p.nc.neq/2; 
    upde=upde(1:p.np*neq); 
else 
    neq=p.nc.neq; 
end
bc=p.fuha.bc(p,p.u); 


