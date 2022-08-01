function r=G(swift_hohenberg_equation,u) 
% rhs for Swift–Hohenberg equation(in the form of 2nd order system)
% also see nonlinearf.m script
% INPUTS
% swift_hohenberg_equation: pde object
%   - nu: dimension of u-components (should euqal 2n)
% u: solution of 2nd order pde system such that u=[u1;u2;par]
%
% OUTPUTs
% r: right hand side of Swift–Hohenberg equation

f=nonlinearf(swift_hohenberg_equation,u); 
r=swift_hohenberg_equation.mat.K*u(1:swift_hohenberg_equation.nu)-swift_hohenberg_equation.mat.M*f; 