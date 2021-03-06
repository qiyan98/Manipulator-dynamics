%% Acceleration tracking results
figure;
ha = tightPlots(1,1,6.85,[5 4],[0.9 0.8],[0.8 0.3],[0.6 0.1],'inch');
cellstr = cell(2,1);
xpos = 0.48; ypos = -0.12;
ha_linewidth = 1.5;
hold on; grid on;

% desired motion
n_plot = 1;
Temp1 = @(t) FunThetaDdot(t);
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
plot(data{m,j,1},data{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
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
plot(data{m,j,1},data{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
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
plot(data{m,j,1},data{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
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
plot(data{m,j,1},data{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of CTC
y_lim = [-4 4.5];
ylim(y_lim);
line([6.8 6.8],y_lim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
line([7.4 7.4],y_lim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);

load SixDOF_simData_ctc.mat;
legend('Desired Motiont','AFNTSM','AFNTSM-SF','LSMC','CTC','location','southwest','numcolumns',2);
% legend('AFNTSM','LSMC','CTC','AFNTSM-SF','location','best','numcolumns',2);
box on;
text(xpos,ypos,'\it (a)','fontweight','bold','units','normalized','fontsize',ha_fontsize);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Acceleration (rad/s^2)','fontsize',ha_fontsize);
a1 = annotation('textarrow',[0.55 0.695],0.92*[1 1],'String','Temperature rise ');
a1.FontSize = ha_fontsize;
a1.Color = 'b';
a2 = annotation('textarrow',[0.55 0.75],0.89*[1 1],'String','External shock ');
a2.FontSize = ha_fontsize;
a2.Color = 'k';

%% Acceleration tracking error results
figure;
% ha = tightPlots(2,2,6.85*0.9,[1 1],[0.6 0.5],[0.7 0.2], [0.7 0.2],'inch');
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
l0 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm_lt = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,3}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM_LT
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Acceleration tracking error (rad/s^2)','fontsize',ha_fontsize);
y_lim = [-1.5 1.6];
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
l0 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_afntsm = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,3}))];
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
% end of AFNTSM-SF
ylim([-3 4]);
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
l0 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_lsmc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,3}))];
n_plot = n_plot + 1;
% end of LSMC
y_lim2 = y_lim;
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Acceleration tracking error (rad/s^2)','fontsize',ha_fontsize);
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
l0 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
max_error_ctc = [max(abs(control_error{m,j,1})) max(abs(control_error{m,j,3}))];
n_plot = n_plot + 1;
% end of CTC
ylim(y_lim2);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('CTC','location','best');
set(gca,'child',[l0 l1 l2]);

axes(ha(3));
text(xpos,ypos,'\it (b)','fontweight','bold','units','normalized','fontsize',ha_fontsize);

%% Disturbance rejection
figure;
ha = tightPlots(2,1,6.85*0.9,[5 2.5],[0.5 0.4],[0.6 0.2], [0.7 0.2],'inch');
%%%%% first plot %%%%%
axes(ha(1));
hold on; box on;grid on;
load SixDOF_simData_afntsm_lt.mat
n_plot = 1;
n_points = 1200;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l00 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_afntsm_sf.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l01 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_lsmc.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l02 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_ctc.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l03 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

xlim([6.795 6.82]);
% xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Acceleration tracking error (rad/s^2)','fontsize',ha_fontsize-2);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
legend('AFNTSM','AFNTSM-SF','LSMC','CTC','location','best','numcolumns',2);
set(gca,'child',[l03 l02 l01 l00 l1]);
a1 = annotation('textarrow',[0.33 0.285],0.75*[1 1],'String','Temperature rise ');
a1.FontSize = ha_fontsize;
a1.Color = 'b';

%%%%% second plot %%%%%
axes(ha(2));
hold on; box on; grid on;
load SixDOF_simData_afntsm_lt.mat
n_plot = 1;
n_points = 100;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l00 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_afntsm_sf.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l01 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_lsmc.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l02 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

load SixDOF_simData_ctc.mat
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l03 = plot(data{m,j,1},control_error{m,j,3},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;
xlim([7.35 7.6]);
ylim([-1.6 1.6]);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM','AFNTSM-SF','LSMC','CTC','location','best','numcolumns',2);
set(gca,'child',[l03 l02 l01 l00 l2]);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Acceleration tracking error (rad/s^2)','fontsize',ha_fontsize-2);
a2 = annotation('textarrow',[0.33 0.285],0.22*[1 1],'String','External shock start');
a2.FontSize = ha_fontsize;
a2.Color = 'k';
a3 = annotation('textarrow',[0.58 0.625],0.15*[1 1],'String','External shock end ');
a3.FontSize = ha_fontsize;
a3.Color = 'k';