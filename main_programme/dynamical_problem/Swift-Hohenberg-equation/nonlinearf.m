function f=nonlinearf(swift_hohenberg_equation,u) 
% SH "nonlinearity" for the 2nd-order system formulation
% f=nonlinearf(swift_hohenberg_equation,u)
%
% INPUTS
% swift_hohenberg_equation: pde object
%   - nu: dimension of u without paramters(should euqal 2n)
% u: solution of 2nd order pde system such that u=[u1;u2;par]
%
% OUTPUTs
% f: nonlinear part of the 2nd order pde system such that f=[f1;f2]

% split parameters from u
par=u(swift_hohenberg_equation.nu+1:end);
lam=par(1); 
nu=par(2);

% split u1 and u2
n=swift_hohenberg_equation.nu/2; 
u1=u(1:n); 
u2=u(n+1:2*n);

% construct the nonlinear part in equation
f1=(lam-1)*u1+nu*u1.^2-u1.^3-2*u2; 
f2=0*u2; 
f=[f1;f2]; 