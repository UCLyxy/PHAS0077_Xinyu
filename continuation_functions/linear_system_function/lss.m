function [x,p]=lss(A,b,p)
% LSS: default LinearSystemSolver, i.e. \ 
%
%  [x,p]=lss(A,b,p)
%

x=A\b; 
