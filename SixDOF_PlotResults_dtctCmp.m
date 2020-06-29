%% Positioning tracking error results
figure;
ha = tightPlots(2,2,6.85*0.9,[1 1],[0.6 0.4],[0.7 0.2], [0.6 0.2],'inch');
xpos = 1.0; ypos = -0.21;
ha_linewidth = 1.5;
n_plot = 1;

%%%%% first subplot %%%%%
axes(ha(1));
hold on; grid on; box on;
% CT-50dB-AFNTSM
load('.mat data backup\CT-50dB\SixDOF_simData_afntsm_lt.mat');
% pos.
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
clear tmpvec tmpindices vec_markerindices;
% end of pos.
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Position tracking error (rad)','fontsize',ha_fontsize);
y_lim_pos = [-2e-3 6e-3];
ylim(y_lim_pos);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM (continuous time)','location','south');
set(gca,'children',[l0 l1 l2]);

%%%%% third plot %%%%%
axes(ha(3)); box on; hold on; grid on;
% vel.
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
n_plot = 3;
l0 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
% end of vel.
y_lim_vel = [-0.04,0.08];
ylim(y_lim_vel);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Velocity tracking error (rad/s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM (continuous time)','location','south');
set(gca,'children',[l0 l1 l2]);


%%%%% second subplot %%%%%
axes(ha(2));
hold on; grid on; box on;
% DT-50dB-AFNTSM
load('.mat data backup\DT-30dB\SixDOF_simData_afntsm_lt.mat');
% pos.
n_points = 10;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
n_plot = 2;
l0 = plot(data{m,j,1},control_error{m,j,1},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
clear tmpvec tmpindices vec_markerindices;
% end of pos.
ylim(y_lim_pos);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM (discrete time)','location','south');
set(gca,'children',[l0 l1 l2]);

%%%%% fourth plot %%%%%
axes(ha(4)); box on; hold on; grid on;
% vel.
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
n_plot = 4;
l0 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
% end of vel.
ylim(y_lim_vel);
xlabel('Time (s)','fontsize',ha_fontsize);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM (discrete time)','location','south');
set(gca,'children',[l0 l1 l2]);

axes(ha(3));
text(xpos,ypos,'\it (a)','fontweight','bold','units','normalized','fontsize',ha_fontsize);
%% plot magnified curve
n_plot = 2;

load('.mat data backup\CT-50dB\SixDOF_simData_afntsm_lt.mat');
annotation('rectangle',[.12 .17 .091 .16],'Color','red','linewidth',ha_linewidth);
annotation('arrow',[.21,.22],[.31,.37],'units','normalized','Color','red');
h3=axes('position',[0.16 0.38 0.2 0.1]);
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
n_plot = n_plot + 1;
axes(h3);
plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
y_lim3 = [-3e-2,3e-2];
xlim([0.5,3]); ylim(y_lim3);
set(gca,'xtick',[0.5:0.5:3]);

load('.mat data backup\DT-30dB\SixDOF_simData_afntsm_lt.mat');
annotation('rectangle',[.59 .17 .091 .16],'Color','red','linewidth',ha_linewidth);
annotation('arrow',[.68,.69],[.32,.37],'units','normalized','Color','red');
h4=axes('position',[0.625 0.38 0.2 0.1]);
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand);
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
n_plot = n_plot + 1;
axes(h4);
plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
xlim([0.5,3]); ylim(y_lim3);
set(gca,'xtick',[0.5:0.5:3]);

%% Disturbance rejection
figure;
ha = tightPlots(2,1,6.85*0.9,[5 2.5],[0.5 0.4],[0.8 0.2], [0.7 0.2],'inch');
%%%%% first plot %%%%%
axes(ha(1));
hold on; box on;grid on;
load('.mat data backup\CT-50dB\SixDOF_simData_afntsm_lt.mat');
n_plot = 1;
n_points = 1200;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l00 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 2;
clear tmpvec tmpindices vec_markerindices;

load('.mat data backup\DT-30dB\SixDOF_simData_afntsm_lt.mat');
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l01 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

xlim([6.795 6.82]);
% xlabel('Time (s)','fontsize',ha_fontsize);
ylim([-2e-2,2e-2]);
ylabel('Velocity tracking error (rad/s)','fontsize',ha_fontsize-2);
l1 = line([6.8 6.8],ylim,'color','blue','linestyle','--','linewidth',0.8*ha_linewidth);
legend('AFNTSM (continuous time)','AFNTSM (discrete time)','location','southwest','numcolumns',2);

set(gca,'child',[l01 l00 l1]);
a1 = annotation('textarrow',[0.33 0.285],0.85*[1 1],'String','Temperature rise ');
a1.FontSize = ha_fontsize;
a1.Color = 'b';

%%%%% second plot %%%%%
axes(ha(2));
hold on; box on; grid on;
load('.mat data backup\CT-50dB\SixDOF_simData_afntsm_lt.mat');
n_plot = 1;
n_points = 100;
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l00 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 2;
clear tmpvec tmpindices vec_markerindices;

load('.mat data backup\DT-30dB\SixDOF_simData_afntsm_lt.mat');
for n_tmp = 1:n_points
    tmpvec = data{m,j,1} - (n_tmp-1+rand)*MaxTime/n_points;
    tmpvec = abs(tmpvec);
    [tmpvec,tmpindices] = sort(tmpvec);
    vec_markerindices(n_tmp) = tmpindices(1);
end
l01 = plot(data{m,j,1},control_error{m,j,2},cell_linespec{n_plot},'color',cell_linecolor{n_plot},'MarkerIndices',vec_markerindices,'linewidth',ha_linewidth);
n_plot = n_plot + 1;
clear tmpvec tmpindices vec_markerindices;

xlim([7.35 7.6]);
ylim([-3e-2,4e-2]);
l2 = line([7.4 7.4],ylim,'color','black','linestyle','-.','linewidth',0.8*ha_linewidth);
legend('AFNTSM (continuous time)','AFNTSM (discrete time)','location','southwest','numcolumns',2);
set(gca,'child',[l01 l00 l2]);
xlabel('Time (s)','fontsize',ha_fontsize);
ylabel('Velocity tracking error (rad/s)','fontsize',ha_fontsize-2);
a2 = annotation('textarrow',[0.33 0.285],0.18*[1 1],'String','External shock start');
a2.FontSize = ha_fontsize;
a2.Color = 'k';
a3 = annotation('textarrow',[0.58 0.625],0.21*[1 1],'String','External shock end ');
a3.FontSize = ha_fontsize;
a3.Color = 'k';

xpos = 0.45; ypos = -0.25;
text(xpos,ypos,'\it (b)','fontweight','bold','units','normalized','fontsize',ha_fontsize);