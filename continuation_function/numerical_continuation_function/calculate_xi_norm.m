function xi_norm=calculate_xi_norm(u,xi,num_equations,weight_aux_eqn)
% calculate_xi_norm: calculate the weighted arclength-continuation xi_norm
% of given vector u
% 
% 
% USAGE:
%       theta_norm=calculate_xi_norm(u,xi,num_equations,xiq)
% INPUT:
% u=given vector(u=[u1;u2;par])
% where u1 is the solution of pde and u2 is the auxiliary variables
% For instance, in swift-hohenberge equation, u2 is the second derivative
% of u1
% xi=value of xi
% num_equations=number of auxiliary equations
% weight_aux_eqn=weight for aux. eqn
% Reference: See dissertation Section 3.3 equation （3.9）
% Called by functions: detect_bifurcation_point.m, detect_fold_point.m
% numerical_continuation.m, initial_step.m, switch_branch.m

length_u=length(u);
% calculate l2 norm of u1
u_l2_norm=norm(u(1:length_u-1-num_equations));
vec_norm=norm(u(length_u-num_equations:length_u)); 
lambda=u(length_u); 
xi_norm=sqrt(xi*u_l2_norm^2+weight_aux_eqn*vec_norm^2+(1-(xi+weight_aux_eqn)/2)*lambda^2);
