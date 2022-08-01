function allen_cahn_equation=oosetfemops(allen_cahn_equation) 
% set stiffness matrix K and mass matrix M
% for Allen-Cahn equation with homogeneous Neumann BCs
[allen_cahn_equation.mat.K,allen_cahn_equation.mat.M,~]=allen_cahn_equation.pdeo.fem.assema(allen_cahn_equation.pdeo.grid,1,1,1); 


