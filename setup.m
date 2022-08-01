%%
% script SETUP: set the path for this Matlab continuation package 
% start using this Matlab program to make numerical continuation
global pphome; pphome=pwd;

fprintf('%s\n',['Matlab continuation program designed by xinyu. Setting library path beginning with ' pphome]); 

addpath(genpath([pphome,'/continuation_functions']));
addpath(genpath([pphome,'/finite_element_geometry']));

format shortG; format compact; 

