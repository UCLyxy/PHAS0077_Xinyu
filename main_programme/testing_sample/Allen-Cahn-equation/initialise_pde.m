function allen_cahn_equation=initialise_pde(lx,nx,par) 
% set up pde object alllen_cahn_equation with some simple geometry, such as:
% domain, boundary condition, grid points for interpolation
%
% INPUTS
% alllen_cahn_equation: initial pde object(structure)
% lx: half range of the domain of solution
% for instance, the domain of solution should be (-lx,lx)
% par: parameters in the equation
%
% OUTPUTS
% alllen_cahn_equation: pde object with some simple geometry

% standard object oriented pde parameters settings
allen_cahn_equation=stanparam;
screenlayout(allen_cahn_equation);

% get the rhs and jacobian of rhs
% save above information into swift_hohenberg_equation
allen_cahn_equation.fuha.sG=@G; 
allen_cahn_equation.fuha.sGjac=@Gu; 

% standard PDE object 1D, yields domain and mesh
pde=stanpdeo1D(lx,2*lx/nx); 
allen_cahn_equation.pdeo=pde;
allen_cahn_equation.np=pde.grid.nPoints; 
allen_cahn_equation.nu=allen_cahn_equation.np;
allen_cahn_equation.sol.xi=1/(allen_cahn_equation.nu); 
% % initial guess, start with the trivial branch u=0 with param appended 
allen_cahn_equation.u=zeros(allen_cahn_equation.np,1); 
allen_cahn_equation.u=[allen_cahn_equation.u; par']; 
allen_cahn_equation.sw.sfem=-1; 
allen_cahn_equation=oosetfemops(allen_cahn_equation); 
% make continuation on par(2)=lambda, up to lammax 
allen_cahn_equation.nc.ilam=2; allen_cahn_equation.nc.lammax=1; 
%  starting and max value for continuation step size
allen_cahn_equation.sol.ds=0.1; allen_cahn_equation.nc.dsmax=0.1; 
allen_cahn_equation.sw.foldcheck=1; % detect and localize folds 
allen_cahn_equation.plot.auxdict={'c','lambda','gamma'}; % parameter names (for axis labels) 