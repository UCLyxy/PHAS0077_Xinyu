function arclength_jac_matrix=gen_ext_mat(p,Gu,G_lambda,tau,xi,xi_aux)
% gen_ext_mat: generate extended (jacobian) matrix for arclength system
% (see dissertation equation(3.13)).
% [Gu,G_lambda;Nu,N_lambda] 
%
% USAGE:
% arclength_jac_matrix=gen_ext_mat(p,Gu,G_lambda,tau,xi,xi_aux)
%
% INPUTS:
% p: object oriented pde object(structure)
% Gu,G_lambda: partial derivatives of the finite-element approximation
% function G(u,lambda) (right hand side of equation)
% tau: tangent vector
% xi: weighting parameter that defines xi-norm
% (see dissertation equation (3.9))
% xi_aux: weigthing of auxiliary equations
% 
% Reference: dissertation equation(3.9, 3.13)
% Called by function: arclength_newton_corrector.m,
% detect_bifurcation_point.m, detect_fold_point.m


if (p.nc.nq>0) % The case that pde has transformed to a higher-order differential system
               % For instance, Swift-Hohenberg equations
               % p.u is made of u1 and u2, u1 is the solution and u2 is the
               % auxiliary varible
    arclength_jac_matrix=[[Gu G_lambda]; [xi*tau(1:p.nu)' xi_aux*tau(p.nu+1:p.nu+p.nc.nq)' (1-(xi+xi_aux))*tau(p.nu+p.nc.nq+1)]];  
else
    arclength_jac_matrix=[[Gu G_lambda]; [xi*tau(1:p.nu)' (1-xi).*tau(p.nu+p.nc.nq+1)]]; 
end

