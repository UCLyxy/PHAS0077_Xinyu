function swift_hohenberg_equation=initialise_pde(swift_hohenberg_equation,nx,lx,par) 
% Swift-Hohenberg as 2nd order system, with 
% singular M=diag(M,0) mass-matrix to compute correct eigenvalues
% set up pde object swift_hohenberg_equation with some simple geometry, such as:
% domain, boundary condition, grid points for interpolation
%
% INPUTS
% swift_hohenberg_equation: initial pde object(structure)
% lx: half range of the domain of solution
% for instance, the domain of solution should be (-lx,lx)
% par: parameters in the equation
%
% OUTPUTS
% swift_hohenberg_equation: pde object with some simple geometry

% standard object oriented pde parameters settings
swift_hohenberg_equation=stanparam(swift_hohenberg_equation);
% 2nd order pde system
swift_hohenberg_equation.nc.neq=2;% 2-component (u1,u2) system
swift_hohenberg_equation.nc.neig=20; 
% get the rhs and jacobian of rhs
% save above information into swift_hohenberg_equation
swift_hohenberg_equation.fuha.sG=@G; 
swift_hohenberg_equation.fuha.sGjac=@Gu; 
% continuation setting
% set up continuation step length
swift_hohenberg_equation.sol.ds=0.01; 
swift_hohenberg_equation.sol.dsmax=0.1; 
swift_hohenberg_equation.pm.resfac=1e-3;
%bifurcation and fold point check
swift_hohenberg_equation.sw.bifcheck=2; 
swift_hohenberg_equation.sw.spcont=0; 
swift_hohenberg_equation.sw.spcalc=1; 
swift_hohenberg_equation.sw.sfem=-1; 
swift_hohenberg_equation.sw.foldcheck=1;
% standard PDE object 1D, yields domain and mesh
pde=stanpdeo1D(lx,2*lx/nx); 
swift_hohenberg_equation.Om=2*lx; 
swift_hohenberg_equation.np=pde.grid.nPoints;  
swift_hohenberg_equation.pdeo=pde;
% pde with a solution u with dimension 2n
% u=[u1;u2;par]
swift_hohenberg_equation.nu=swift_hohenberg_equation.np*swift_hohenberg_equation.nc.neq; 
swift_hohenberg_equation.sol.xi=1/swift_hohenberg_equation.nu;
swift_hohenberg_equation.nc.lammin=-4; 
swift_hohenberg_equation.nc.lammax=2; 
swift_hohenberg_equation.nc.ds=0.01; 
swift_hohenberg_equation.nc.dsmax=0.1;
% initial guess, start with the trivial branch u=0
u=0*ones(swift_hohenberg_equation.np,1); 
v=u; 
u0=[u v]; 
swift_hohenberg_equation.u=u0(:); 
swift_hohenberg_equation.u=[swift_hohenberg_equation.u; par]; 
swift_hohenberg_equation.nc.ilam=1; 
swift_hohenberg_equation.file.smod=5; 
swift_hohenberg_equation=setfemops(swift_hohenberg_equation);
swift_hohenberg_equation.nc.nsteps=100;