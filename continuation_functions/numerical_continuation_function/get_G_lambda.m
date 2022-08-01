function G_lambda=get_G_lambda(p,u,r)
% get_G_lambda: get the finite difference approximation of derivative of G
%  with respect to primary parameter defined in p.nc.ilam
%
%  USAGE:
%  G_lambda=get_G_lambda(p,u,r)
%
%  INPUTS:
%   *p: object oriented pde object
%   *u: solution vector of pde system(u=[u1;u2;par])
%   *r: residual computed by residual_pde.m 
%
%  RETURNS:
%   *G_lambda: finite difference approximation of derivative of G
%              w.r.t parameter lambda
%
% Reference: See dissertaion: total derivative equation (3.3)

u(p.nu+p.nc.ilam(1))=u(p.nu+p.nc.ilam(1))+p.nc.del; % increment primary bif. param.
r1=calculate_residual(p,u); % size(r1), size(r)
G_lambda=(r1-r)/p.nc.del; 