%% command line interface for SH equation as consistent 2nd order system 
% Please run this script cell by cell to see what happens at each step
keep pphome; 
close all; 
swift_hohenberg_equation=[]; 
%% initialise and set up zero-branch (trivial solution) on long domain
% set domain (-10pi,10pi)
lx=10*pi; 
nx=round(50*lx);

% set parameters
lam=-0.05; 
nu=2; 
par=[lam; nu];

swift_hohenberg_equation=initialise_pde(swift_hohenberg_equation,nx,lx,par); 
huclean(swift_hohenberg_equation);
% set name of the trivial branch
swift_hohenberg_equation=set_file_name(swift_hohenberg_equation,'Example2/trivial');
% bifurcation detect
swift_hohenberg_equation=find_bifurcation_point(swift_hohenberg_equation,4);
% All these four bifurcation points (bpt) would be recorded in the
% directory "trivial"
%% first Turing-branch
% switch branch at the first detected bifurcation point
% and then make continuation
swift_hohenberg_equation=switch_branch('Example2/trivial','bpt1','Example2/branch1',0.01); 
swift_hohenberg_equation=numerical_continuation(swift_hohenberg_equation,40);
% numerical_continuation itself contains the function of detecting
% bifurcation points along the proceeding of continuation
%% bifurcation on bifurcation 
swift_hohenberg_equation=switch_branch('Example2/branch1','bpt1','Example2/subranch1',0.01); 
% note that here 'bpt1' is the first bifurcation point detected on
% 'branch1' rather than 'trivial'
swift_hohenberg_equation.nc.dsmax=0.1; 
swift_hohenberg_equation=numerical_continuation(swift_hohenberg_equation,180);
%% second sub branch
swift_hohenberg_equation=switch_branch('Example2/branch1','bpt2','Example2/subranch2',-0.01); 
swift_hohenberg_equation.nc.dsmax=0.05; 
swift_hohenberg_equation.nc.tol=1e-6; 
swift_hohenberg_equation.sw.bifcheck=0; 
swift_hohenberg_equation.pm.resfac=1e-2; 
swift_hohenberg_equation=parallel_multi_continuation(swift_hohenberg_equation,250); 
%% Solution plotting on branch1
plotsol('Example2/branch1','pt5');
plotsol('Example2/branch1','pt30');
plotsol('Example2/branch1','pt40');
%% Solution plotting on subranch1
plotsol('Example2/subranch1','pt5');
plotsol('Example2/subranch1','pt90');
plotsol('Example2/subranch1','pt180');
%% Solution plotting on subranch2
plotsol('Example2/subranch2','pt5');
plotsol('Example2/subranch2','pt100');
plotsol('Example2/subranch2','pt250');
