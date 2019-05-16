
bayes_root = pwd; 

cd '../../grasp_m_barebones'

grasp_ini;
grasp_startup;
disp('Adding override paths for Bayesian stuff');
%

addpath(fullfile(bayes_root,'GRASP_overrides_v8.07')); 
addpath(genpath(fullfile(bayes_root,'minFunc_2012')));
addpath(bayes_root);
disp('Run grasp.m as usual now');