%% create a new plot for control input - separated
figure;
ha = tightPlots(2,2,6.85*0.9,[1 1],[0.6 0.5],[0.7 0.2], [0.7 0.2],'inch');
xpos = 1.0; ypos = -0.21;
ha_linewidth = 1.5;
n_plot = 1;

%%%%% first subplot %%%%%
axes(ha(1));
hold on; grid on; box on;
% AFNTSM
load SixDOF_simData_afntsm_lt.mat;
n_points = 80;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,5},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm_lt = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM_LT
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Control input (Nm)','fontsize',ha_fontsize);
y_lim = [-0 6e3];
x_lim = [7.3 7.8];
xlim(x_lim); ylim(y_lim);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM','location','best');
set(gca,'children',[l0 l1 l2]);

%%%%% second subplot %%%%%
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
l0 = plot(data{m,j,1},data{m,j,5},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM-SF
xlim(x_lim); ylim(y_lim);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM-SF','location','best');
set(gca,'children',[l0 l1 l2]);

%%%%% third plot %%%%%
axes(ha(3)); box on; hold on; grid on;
% LSMC
load SixDOF_simData_lsmc;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,5},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_lsmc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of LSMC
y_lim2 = y_lim;
% y_lim2 = [-7e-3 0.025];
xlim(x_lim);
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Control input (Nm)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('LSMC','location','best');
set(gca,'children',[l0 l1 l2]);

%%%%% fourth plot %%%%%
axes(ha(4)); box on; hold on; grid on;
% CTC
load SixDOF_simData_ctc;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},data{m,j,5},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_ctc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of CTC
xlim(x_lim);
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('CTC','location','best');
set(gca,'children',[l0 l1 l2]);

axes(ha(3));
% text(xpos,ypos,'\it (b)','fontweight','bold','units','normalized','fontsize',ha_fontsize);