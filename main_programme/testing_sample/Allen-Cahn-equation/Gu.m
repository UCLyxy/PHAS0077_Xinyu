function Gu=Gu(p,u)
% Jacobian of the rhs for Allen-Cahn equation
% with homogeneous Neumann BCs
% also see G.m script

% split u into parameters and PDE variables 
par=u(p.nu+1:end); u=u(1:p.nu);
% local derivative of 'nonlinearity'  
fu=par(2)+3*u.^2-5*par(3)*u.^4; 
% put derivatives into (sparse) matrix 
Fu=spdiags(fu,0,p.nu,p.nu);  
% the Jacobian matrix 
Gu=par(1)*p.mat.K-p.mat.M*Fu;    