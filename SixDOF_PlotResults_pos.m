%% Position tracking results
figure;
ha = tightPlots(1,1,6.85,[5 4],[0.9 0.8],[0.8 0.3],[0.6 0.1],'inch');
cellstr = cell(2,1);
xpos = 0.48; ypos = -0.12;
ha_linewidth = 1.5;
hold on; grid on;

% desired motion
n_plot = 1;
Temp1 = @(t) FunTheta(t);
fplot(Temp1,[0 MaxTime],'r--','linewidth',ha_linewidth);
% end of desired motion

% AFNTSM
load SixDOF_simData_afntsm_lt.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
plot(data{m,j,1},data{m,j,2}(:,1),cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM

% AFNTSM-SF
load SixDOF_simData_afntsm_sf.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
plot(data{m,j,1},data{m,j,2}(:,1),cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM-SF

% LSMC
load SixDOF_simData_lsmc.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
plot(data{m,j,1},data{m,j,2}(:,1),cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of LSMC

% CTC
load SixDOF_simData_ctc.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
plot(data{m,j,1},data{m,j,2}(:,1),cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of CTC
y_lim = [1.05 1.4];
line([6.8 6.8],y_lim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
line([7.4 7.4],y_lim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);

load SixDOF_simData_ctc.mat;
legend('Desired Motiont','AFNTSM','AFNTSM-SF','LSMC','CTC','location','southwest','numcolumns',2);
% legend('AFNTSM','LSMC','CTC','AFNTSM-SF','location','best','numcolumns',2);
box on;
text(xpos,ypos,'\it (a)','fontweight','bold','units','normalized','fontsize',ha_fontsize);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Position (rad)','fontsize',ha_fontsize);
a1 = annotation('textarrow',[0.55 0.695],0.92*[1 1],'String','Temperature rise ');
a1.FontSize = ha_fontsize;
a1.Color = 'b';
a2 = annotation('textarrow',[0.55 0.75],0.89*[1 1],'String','External shock ');
a2.FontSize = ha_fontsize;
a2.Color = 'k';
%% Positioning tracking error results
figure;
ha = tightPlots(2,2,6.85*0.9,[1 1],[0.6 0.4],[0.7 0.2], [0.6 0.2],'inch');
xpos = 1.0; ypos = -0.21;
ha_linewidth = 1.5;
n_plot = 1;

%%%%% first subplot %%%%%
axes(ha(1));
hold on; grid on; box on;
% AFNTSM
load SixDOF_simData_afntsm_lt.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm_lt = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM_LT
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Position tracking error (rad)','fontsize',ha_fontsize);
y_lim = [-4e-3 8e-3];
ylim(y_lim);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM','location','best');
set(gca,'children',[l0 l1 l2]);

%%%%% second subplot %%%%%
axes(ha(2));
hold on; grid on; box on;
% AFNTSM-SF
load SixDOF_simData_afntsm_sf.mat;
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM-SF
ylim(y_lim);
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
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_lsmc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of LSMC
y_lim2 = y_lim;
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Position tracking error (rad)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('LSMC','location','best');
set(gca,'children',[l0 l1 l2]);

%%%%% fourth plot %%%%%
axes(ha(4)); box on; hold on; grid on;
% CTC
load SixDOF_simData_ctc;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_ctc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,2}))];
n_plot = n_plot + 1;
% end of CTC
% y_lim2 = [-6e-3,9e-3];
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('CTC','location','best');
set(gca,'children',[l0 l1 l2]);

axes(ha(3));
text(xpos,ypos,'\it (b)','fontweight','bold','units','normalized','fontsize',ha_fontsize);