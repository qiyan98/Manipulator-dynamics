close all;
clear;
tic;
addpath('dynamics_derivation');
addpath('visualization_tools');
addpath('visualization_tools\export_fig');

syms t real;
MaxTime = 10;
FunTheta = 70/57.3 + 4*(1+(exp(t)/(exp(t)+30)))*cos(3*t*(1+0.5*(exp(t)/(exp(t)+30))))/57.3; % Trajectory is set here!!!

this_tem = 45; this_t3 = 0;  
save SixDOF_traj.mat FunTheta MaxTime this_tem this_t3;

% SixDOF_TrajSim1(FunTheta,MaxTime); % Computed torque without friction
% SixDOF_simulation([40 50 60],[0]./57.3,FunTheta,MaxTime,'openloop','true');

% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'pid','error');
% 
% clear; load SixDOF_traj.mat;
% disp('Simulation on CTC...');
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'ctc','error');
% 
% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'fntsm','error');
% 
% 
% clear; load SixDOF_traj.mat;
% disp('Simulation on LSMC...');
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'lsmc','error');
% % 
% clear; load SixDOF_traj.mat;
% disp('Simulation on AFNTSM-SF...');
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'afntsm_sf','error');
% % % 
clear; load SixDOF_traj.mat;
disp('Simulation on AFNTSM-LT...');
SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'afntsm_lt','error');

% 
% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'hinf','error');
% SixDOF_PlotResults;

SixDOF_PlotResults;
toc;

% validation
% SixDOF_simulation_validation(FunTheta,MaxTime);

% SixDOF_FricPlot;