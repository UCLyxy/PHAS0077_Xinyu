function determainant=get_determainant(p)
% get_determainant: get determinant of the jacobian of the
% extended system(See gen_ext_mat.m for details)
% 
% See dissertation equation (3.13)
%
% USAGE:
%  determainant=get_determainant(p)
%
% In case p.nc.neigdet=0 use LU-decomposition.
% 



r=calculate_residual(p,p.u);
[Gu,G_lambda]=get_info_derivative(p,p.u,r);
arclength_jac_matrix=gen_ext_mat(p,Gu,G_lambda,p.tau,p.sol.xi,p.sol.xiq);

if p.nc.neigdet==0
   try; [L,U,P,Q]=lu(arclength_jac_matrix);
   determainant=full(prod(sign(diag(U)))*det(P)*det(Q));
   catch; [L,U]=lu(arclength_jac_matrix);
   determainant=full(prod(sign(diag(U)))); %*det(P)*det(Q));
   end 
else
  if p.sw.eigsstart==1; vs=size(arclength_jac_matrix,1); p.sw.evopts.v0=ones(vs,1)/vs; end 
  [V,mu]=eigs(arclength_jac_matrix,p.nc.neigdet,p.nc.eigref,p.sw.evopts); muv=mu*ones(1,p.nc.neigdet)'; 
  determainant=prod(sign(real(muv))); 
end
