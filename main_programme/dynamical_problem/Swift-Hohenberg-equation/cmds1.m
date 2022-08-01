%% command line interface for SH equation as consistent 2nd order system 
% Please run this script cell by cell to see what happens at each step
keep pphome; 
close all; 
swift_hohenberg_equation=[]; 
%% initialise and set up zero-branch (trivial solution) on long domain
% set domain (-10pi,10pi)
lx=20*pi; 
nx=round(30*lx);
% set parameters
lam=-0.05;
nu=2; 
par=[lam;nu]; 

swift_hohenberg_equation=initialise_pde(swift_hohenberg_equation,nx,lx,par);
swift_hohenberg_equation.plot.pmod=10;
% set branch name
swift_hohenberg_equation=set_file_name(swift_hohenberg_equation,'Example1/trivial'); 
% detect bifurcation points
swift_hohenberg_equation=find_bifurcation_point(swift_hohenberg_equation,4);
% All these four bifurcation points (bpt) would be recorded in the
% directory "trivial"
%% Turing-branches 
swift_hohenberg_equation=switch_branch('Example1/trivial','bpt1','Example1/branch1',0.01); 
title("bifurcation diagram");
swift_hohenberg_equation.sw.bifcheck=0; 
swift_hohenberg_equation=numerical_continuation(swift_hohenberg_equation,40);