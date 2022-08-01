function r=G(allen_cahn_equation,u)
% rhs for Allen-Cahn equation with homogeneous Neumann BCs
% INPUTS
% allen_cahn_equation: pde object
% u: solution of pde such that u=[u;par]
%
% OUTPUTs
% r: right hand side of Allen-Cahn equation

par=u(allen_cahn_equation.nu+1:end); 
% split u into parameters and PDE variables 
u=u(1:allen_cahn_equation.nu);
% "nonlinearity", i.e., everything but diffusion 
f=par(2)*u+u.^3-par(3)*u.^5; 
% the residual (rhs)
r=par(1)*allen_cahn_equation.mat.K*u-allen_cahn_equation.mat.M*f; 