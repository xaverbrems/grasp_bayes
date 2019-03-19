
bayes_root = pwd; 

cd '../../Grasp_v814/grasp_m_barebones'

grasp_startup;
disp('Adding override paths for Bayesian stuff');
%

addpath(fullfile(bayes_root,'GRASP_overrides_v8.07')); 
addpath(genpath(fullfile(bayes_root,'minFunc_2012'))); 

grasp