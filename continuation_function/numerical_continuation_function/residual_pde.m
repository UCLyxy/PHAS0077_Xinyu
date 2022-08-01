function r=residual_pde(p,u)
% residual_pde: return pde-part of residual
% which is just right hand side of our finite_element approximation
% function G(u,u(lambda))
%  USAGE:
%   r=residual_pde(p,u)
%
%  INPUTS:
%   *p: object oriented pde object
%   *u: solution vector of pde system(u=[u1;u2;par])
%
%  RETURNS:
%   *r: residual for PDE-part
%
%  Called by function: calculate_residual

r=p.fuha.sG(p,u);