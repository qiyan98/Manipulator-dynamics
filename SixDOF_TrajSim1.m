function SixDOF_TrajSim1(FunTheta,MaxTime)
% Numerical Dynamic simulation - inverse dynamics/computed torque without friction

% preparation for equation solver
SixDOF_dynamics_simplified; % load necessary dynamic variables
SixDOF_FunGeneration; % load necessary anonymous functions

% Input theta_3
% t3 = 0;
% FunDesTorque2 = @(t)FunDesTorque2(t,t3);
% FunG2 = @(t2) FunG2(t2,t3);

% Plot results
figure;
subplot(3,1,1);
% hold on;
% fplot(FunDesTorque2,[0 MaxTime]);
% line(xlim,[-6310.98 -6310.98],'linestyle','--','color','red');
% line(xlim,[6310.98 6310.98],'linestyle','--','color','red');
% xlabel('Time (s)');
% ylabel('Motor torque (Nm)');
% legend('Torque by motor 2');
hold on;
fplot(FunM22,[-pi 70/57.3]);
% line(xlim,[-6310.98 -6310.98],'linestyle','--','color','red');
% line(xlim,[6310.98 6310.98],'linestyle','--','color','red');
% xlim([0 2*pi]);
% xticks(0:pi/2:2*pi);
% xticklabels({'0','0.5\pi','\pi','1.5\pi','2\pi'});
xlim([-pi 70/57.3]);
ylim([min(ylim) 550]);
xlabel('\theta_3 (rad)');
ylabel('M_{22} (kg\cdotm^2)');
box on;
% legend('M_{22}');


[t2,t3] = meshgrid(-65/57.3:0.01:140/57.3,-3:0.1:1);
G_mat = zeros(size(t2,1),size(t3,2));
for i = 1:size(t2,1)
    for j = 1:size(t3,2)
        G_mat(i,j) = FunG2(t2(i,j),t3(i,j));
    end
end

subplot(3,1,2:3);
box on;
% [b,h] = contour(t2,G_mat,t3,'levelstep',1);
[b,h] = contour(t2,G_mat,t3,3);
clabel(b,h,'labelspacing',160,'fontsize',8);
c = colorbar;
c.Label.String = '\theta_3 (rad)';
colormap('jet');

% xlim([0 2*pi]);
% xticks(0:pi/2:2*pi);
% xticklabels({'0','0.5\pi','\pi','1.5\pi','2\pi'});
xlabel('\theta_2 (rad)');
ylabel('G_2 (Nm)');
zlabel('\theta_3 (rad)');
set(gcf,'color','white');
% export_fig effect_configuration.eps;
% tightfigadv;

% hold on;
% line(xlim,[-6310.98 -6310.98],'linestyle','--','color','red');
% line(xlim,[6310.98 6310.98],'linestyle','--','color','red');
% xlabel('Position of joint 2 (deg)');
% ylabel('Gravitational Torque \tau_l (Nm)');
% legend('Gravitational torque');

% figure;
% subplot(3,1,1);
% Temp1 = @(t) FunTheta(t).*57.3;
% fplot(Temp1,[0 MaxTime]);
% hold on;
% % line(xlim,[-65 -65],'linestyle','--','color','red');
% % line(xlim,[140 140],'linestyle','--','color','red');
% xlabel('Time/s');
% ylabel('Deg');
% title('Joint 2 Position');
% 
% subplot(3,1,2);
% Temp1 = @(t) FunThetaDot(t).*57.3;
% fplot(Temp1,[0 MaxTime]);
% hold on;
% % line(xlim,[-90 -90],'linestyle','--','color','red');
% % line(xlim,[90 90],'linestyle','--','color','red');
% xlabel('Time/s');
% ylabel('Deg/s');
% title('Joint 2 Velocity');
% 
% subplot(3,1,3);
% Temp1 = @(t) FunThetaDdot(t).*57.3;
% fplot(Temp1,[0 MaxTime]);
% xlabel('Time/s');
% ylabel('Deg/s^2');
% title('Joint 2 Acceleration');
end