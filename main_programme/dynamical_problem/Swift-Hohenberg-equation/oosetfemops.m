function swift_hohenberg_equation=oosetfemops(swift_hohenberg_equation) 
% set stiffness matrix K and mass matrix M
% for Swift–Hohenberg equation as 2nd order system


% scalar stiffness and mass matrix
[K,M,~]=swift_hohenberg_equation.pdeo.fem.assema(swift_hohenberg_equation.pdeo.grid,1,1,1); 

% 2nd order system stiffness 
swift_hohenberg_equation.mat.K=[[0*K -K];[K M]];
% 2nd order system mass matrix (here is singular)
swift_hohenberg_equation.mat.M=[[M 0*M];[0*M 0*M]];  