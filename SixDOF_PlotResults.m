%% Initialization
clear;
rng('shuffle');
% close all;
load SixDOF_simData_ctc.mat;
if mode_ctrl == -1
    return;
end
if strcmp(options,'true')
    plotTrue = 1;
elseif strcmp(options,'error')
    plotTrue = 0;
end

if n_t3 == 1
	cellstr{1,1} = ['Desired Motion'];
    tmp = 2;
    for m=n_t3
        for j=n_tem
            cellstr{tmp,1} = ['T = ',num2str(CurT(j)),'¡æ, ','\theta_3 = ',num2str(Configuration(m)*57.3),'¡ã'];
            tmp = tmp + 1;
        end
    end
elseif n_t3 > 1
    cellstr{1} = ['Desired Motion'];
    tmp = 2;
    for m=1:n_t3
        for j=1:n_tem
            cellstr{tmp,1} = ['T = ',num2str(CurT(j)),'¡æ, \theta_3 = ',num2str(Configuration(m)*57.3),'¡ã'];
            tmp = tmp + 1;
        end
    end
end
if plotTrue == 0 % error mode
    cellstr_origin = cellstr;
    cellstr_1 = cellstr(1);
    cellstr = cellstr(2:end);
end
cell_linespec = cell(10,1);
% cell_linespec{1} = '-o'; cell_linespec{2} = '-*';
% cell_linespec{3} = '-d'; cell_linespec{4} = '-^'; 
cell_linespec{1} = '-o'; cell_linespec{2} = '--*';
cell_linespec{3} = ':d'; cell_linespec{4} = '-.^'; 
cell_linespec{5} = '-.>'; cell_linespec{6} = '-.h';
cell_linespec{7} = '-.';
cell_linecolor = cell(10,1);
cell_linecolor{1} = [0, 0.4470, 0.7410]; cell_linecolor{2} = [0.4660, 0.6740, 0.1880];
cell_linecolor{3} = [0.9290, 0.6940, 0.1250]; cell_linecolor{4} = [0.4940, 0.1840, 0.5560];
cell_linecolor{5} = [0.4660, 0.6740, 0.1880];
n_points = 10;
ha_fontsize = 13;
vec_markerindices = zeros(2,1);
%% Pos
SixDOF_PlotResults_pos();
%% Vel
SixDOF_PlotResults_vel();
%% Acc
% SixDOF_PlotResults_acc();
%% Ctrl_input
% SixDOF_PlotResults_ctrl_input();
%% FricTorq
SixDOF_PlotResults_fricTorq();
%% Estimated parameters
% load SixDOF_simData_afntsm.mat
% UB_est = zeros(10,1);
% UB_real = zeros(10,1);
% for i = 1:length(X_ode)
%     UB_est(i) = X_ode(i,3) + X_ode(i,4)*abs(X_ode(i,1)) + X_ode(i,5)*abs(X_ode(i,2));
%     UB_real(i) = data{m,j,6}(i) + FunD(t_ode(i));
% end
% figure; box on; hold on;
% plot(t_ode,UB_est);
% plot(t_ode,X_ode(:,3));
% plot(t_ode,X_ode(:,4));
% plot(t_ode,X_ode(:,5));
% legend('Estimated UB','a_0','a_1','a_2','location','best');