%% Contribution to friction torque change of various factors
figure;
ha = tightPlots(1,1,6.85,[5 4],[0.9 0.8],[0.8 0.3],[0.7 0.1],'inch');
cellstr = cell(2,1);
xpos = 0.48; ypos = -0.12;
ha_linewidth = 1.5;
hold on; grid on;
n_plot = 1
% AFNTSM
load SixDOF_simData_afntsm_lt.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
plot(data{m,j,1},data{m,j,4},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
plot(data{m,j,1},data{m,j,6},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
% plot(data{m,j,1},abs(data{m,j,7})./(abs(data{m,j,7})+abs(data{m,j,8})),cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
plot(data{m,j,1},data{m,j,7},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
plot(data{m,j,1},data{m,j,8},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM
ylim([-1.2e3,1.5e3]);
line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);

legend('\tau_f','\Delta\tau_f',...
    '\Delta\tau_{fT}','\Delta\tau_{fL}','location','southwest','numcolumns',2);
box on;
% text(xpos,ypos,'\it (a)','fontweight','bold','units','normalized','fontsize',ha_fontsize);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Friction torque (Nm)','fontsize',ha_fontsize);
a1 = annotation('textarrow',[0.55 0.695],0.92*[1 1],'String','Temperature rise ');
a1.FontSize = ha_fontsize;
a1.Color = 'b';
a2 = annotation('textarrow',[0.55 0.75],0.89*[1 1],'String','External shock ');
a2.FontSize = ha_fontsize;
a2.Color = 'k';
%% create a new plot for control input - separated
figure;
ha = tightPlots(2,2,6.85*0.9,[1 1],[0.6 1.0],[0.7 0.2], [0.8 0.2],'inch');
set(gcf,'color','w');
xpos = 1.0; ypos = -0.21;
ha_linewidth = 1.5;
n_plot = 1;

%%%% first subplot %%%%%
axes(ha(1));
hold on; grid on; box on;
% AFNTSM
load SixDOF_simData_afntsm_lt.mat;
n_points = 5;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,4},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm_lt = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM_LT
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('$$\tau_f$$ (Nm)','interpreter','latex','fontsize',ha_fontsize);
y_lim = [-1.5e3 1.5e3];
x_lim = [0 10];
xlim(x_lim); ylim(y_lim);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
set(gca,'children',[l0 l1 l2]);

%%%% second subplot %%%%%
axes(ha(2));
hold on; grid on; box on;
% AFNTSM-SF
load SixDOF_simData_afntsm_sf.mat;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,6},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM-SF
y_lim = [-1500,1500];
xlim(x_lim); ylim(y_lim);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('$$\Delta\tau_f = (\tau_f - \tau_f^*)$$ (Nm)','interpreter','latex','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
set(gca,'children',[l0 l1 l2]);
set(gca,'Xtick',0:2:10)

%%%% third plot %%%%%
axes(ha(3)); box on; hold on; grid on;
% LSMC
load SixDOF_simData_lsmc;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,7},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_lsmc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of LSMC
y_lim2 = [-400,800];
xlim(x_lim);
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel({'$$\Delta\tau_{fT}$$ friction variation','due to temperature (Nm)'},...
    'interpreter','latex','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
set(gca,'children',[l0 l1 l2]);
set(gca,'Xtick',0:2:10)

%%%% fourth plot %%%%%
axes(ha(4)); box on; hold on; grid on;
% CTC
load SixDOF_simData_ctc;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,8},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_ctc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of CTC
xlim(x_lim);
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel({'$$\Delta\tau_{fL}$$ friction variation','due to load (Nm)'},'interpreter','latex',...
    'fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
set(gca,'children',[l0 l1 l2]);
set(gca,'Xtick',0:2:10)

axes(ha(3));
% text(xpos,ypos,'\it (b)','fontweight','bold','units','normalized','fontsize',ha_fontsize);