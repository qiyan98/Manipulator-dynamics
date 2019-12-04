close all;
clear;
tic;
addpath('dynamics_derivation');
addpath('visualization_tools');
addpath('visualization_tools\export_fig');

syms t real;
% FunTheta = 70/57.3 + 0.8*sin(t) + 0.1*sin(2*t); % Trajectory is set here!!!
% FunTheta = 70/57.3 + (6*sin(5*t) + 8*cos(4*t))/57.3; % Trajectory is set here!!!
MaxTime = 10;
% FunTheta = 70/57.3 + 4*(1+(exp(t)/(exp(t)+30)))*cos(3*t*(1+0.5*(exp(t)/(exp(t)+30))))/57.3; % Trajectory is set here!!!
FunTheta = 70/57.3 + 4*(1+(exp(t)/(exp(t)+30)))*cos(3*t*(1+0.5*(exp(t)/(exp(t)+30))))/57.3; % Trajectory is set here!!!

% FunTheta = (100 -  20*(1/(exp(100*t))))/57.3; % Trajectory is set here!!!
% FunTheta = 70/57.3 + 0*t;
% FunTheta = 70/57.3+5*sin(t/100)/57.3;
this_tem = 45; this_t3 = 0;  
save SixDOF_traj.mat FunTheta MaxTime this_tem this_t3;
% FunTheta = 70/57.3 + 5*cos(4*t*(0.2))/57.3; % Trajectory is set here!!!

% SixDOF_TrajSim1(FunTheta,MaxTime); % Computed torque without friction
% SixDOF_simulation([40 50 60],[0]./57.3,FunTheta,MaxTime,'openloop','true');

% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'pid','error');
% 
clear; load SixDOF_traj.mat;
SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'ctc','error');
% 
% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'fntsm','error');
% 
clear; load SixDOF_traj.mat;
SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'lsmc','error');
% % 
clear; load SixDOF_traj.mat;
SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'afntsm_sf','error');
% % 
clear; load SixDOF_traj.mat;
SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'afntsm_lt','error');
% 
% clear; load SixDOF_traj.mat;
% SixDOF_simulation([this_tem],[this_t3]./57.3,FunTheta,MaxTime,'hinf','error');

% SixDOF_PlotResults;
SixDOF_PlotResultsNew;
toc;
% validation
% SixDOF_simulation_validation(FunTheta,MaxTime);

% SixDOF_FricPlot;