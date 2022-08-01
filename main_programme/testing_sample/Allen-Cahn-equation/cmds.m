%% C1, preparations, close windows, clear workspace
% this is a test for the validation of program
% Reference: dissertation Section 4.4
close all; 
keep pphome; 
%% C2: initialize(generic), then specific settings (could also be set in initialise_pde.m)
par=[1 -0.2 1]; 
lx=5; 
nx=30; 
allen_cahn_equaion=initialise_pde(lx,nx,par); 
% set file name of output directory for trivial branch
allen_cahn_equaion=set_file_name(allen_cahn_equaion,'trivial'); 
%% C3: find bifpoints of the trivial branch
allen_cahn_equaion=find_bifurcation_point(allen_cahn_equaion,4);
%% C4: switch to bifurcating branches and continue
allen_cahn_equaion=switch_branch('trivial','bpt1','branch1',0.1); 
allen_cahn_equaion=numerical_continuation(allen_cahn_equaion); 
%% C5: solution plots 
plotsol('branch1','pt5'); 
axis([-5 5 0 1]); 
title("5-th point of constant branch");
ylabel("u(x)");
%% C6: plotting of the exact solution
a1=sqrt(10);
f1=@(t) a1.*sqrt(0.5*(1-sqrt(1+4.*t)));
f2=@(t) a1.*sqrt(0.5*(1+sqrt(1+4.*t)));
t1=linspace(-0.25,0,201);
figure,plot(t1,f1(t1));
xlabel("\lambda");
ylabel("||u||_2");
hold on;
t2=linspace(-0.25,0.3,201);
plot(t2,f2(t2));
title("exact solution of constant branch");
xlim([-0.5 1]);
