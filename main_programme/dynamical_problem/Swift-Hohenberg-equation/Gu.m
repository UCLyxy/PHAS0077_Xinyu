function Gu=Gu(swift_hohenberg_equation,u)
% Jacobian of the rhs for Swift–Hohenberg equation
% (in the form of 2nd order system)
% also see nonlinearf.m script and G.m script
% INPUTS
% swift_hohenberg_equation: pde object
%   - nu: dimension of u-components (should euqal 2n)
% u: solution of 2nd order pde system such that u=[u1;u2;par]
%
% OUTPUTs
% Gu: Jacobian of the right hand side of Swift–Hohenberg equation as
% 2nd order system
[f1u,f1v,f2u,f2v]=njac(swift_hohenberg_equation,u); n=swift_hohenberg_equation.nu/2;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=swift_hohenberg_equation.mat.K-swift_hohenberg_equation.mat.M*Fu; 
end 

function [f1u,f1v,f2u,f2v]=njac(swift_hohenberg_equation,u) 
n=swift_hohenberg_equation.nu/2; u1=u(1:n); par=u(swift_hohenberg_equation.nu+1:end); lam=par(1); nup=par(2); ov=ones(n,1);  
f1u=-(1-lam)*ov+2*nup*u1-3*u1.^2; f1v=-2*ov;
f2u=0*ov; f2v=0*ov; 
end