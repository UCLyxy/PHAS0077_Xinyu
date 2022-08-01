function [Gu,G_lambda]=get_info_derivative(p,u,r)
% get_info_derivative: return jacobian Gu and derivative Glam depending on p.sw.jac 
%
%  USAGE:
%  [Gu,G_lambda]=get_info_derivative(p,u,r)
%
%  INPUTS:
%   *p: object oriented pde object
%   *u: solution vector of pde system(u=[u1;u2;par])
%   *r: residual computed by residual_pde.m
%
%  RETURNS:
%   *Gu: partial derivatives of the finite-element approximation function G(u,lambda) 
%        (right hand side of equation) with respect to u
%   *G_lambda: finite difference approximation of derivative of G
%              w.r.t parameter lambda
%  
%  Reference: See dissertaion: total derivative equation (3.3)

Gu=get_Gu(p,u,r);
G_lambda=get_G_lambda(p,u,r); 
