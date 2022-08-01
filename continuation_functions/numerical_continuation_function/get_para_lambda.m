function lambda=get_para_lambda(p,varargin)
% GET_PARA_LAMBDA: return value of current primary parameter.
%
%  Important settings: p.nc.ilam
% 
% Called by function: detect_bifurcation_point.m, detect_fold_point,
%                      numerical_continuation.m

if ~isempty(varargin) 
    u=varargin{1}; 
else 
    u=p.u; 
end
    lambda=u(p.nu+p.nc.ilam(1));
end
