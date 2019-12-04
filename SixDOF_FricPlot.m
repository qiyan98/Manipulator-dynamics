% Plot graphs for the complex friction model with respect to
% temperature and load torque

% % FunFriction(tau_a,theta_dot,tau_load,T)
close all;
addpath('visualization_tools');
addpath('visualization_tools\export_fig');
%%
[theta_dot,tau_load] = meshgrid(0:1e-2:1,0:1000:6000);
Fric = zeros(size(theta_dot,1),size(tau_load,2));
for n_x = 1:size(theta_dot,1)
    for n_y = 1:size(tau_load,2)
        Fric(n_x,n_y) = FunFriction(tau_load(n_x,n_y),theta_dot(n_x,n_y),40);
    end
end

figure;
% [b,h] = contour(theta_dot,Fric,tau_load,'showtext','on','LabelSpacing',500);
[b,h] = contour(theta_dot,Fric,tau_load,4);
clabel(b,h,'labelspacing',500,'fontsize',7);
c = colorbar;
c.Label.String = 'Load torque (Nm)';
colormap('winter');
ylim([200 1200]);
x = [0.4 0.345];
y = [0.23 0.39];
annotation('textarrow',x,y,'String','Load torque increases','linewidth',1);
xlabel('$$\dot{\theta}(rad/s)$$','interpreter','latex');
ylabel('$$\tau_f(Nm)$$','interpreter','latex');
set(gcf,'color','white');
% export_fig Fric_Load_curve.eps;
% title('Friction curves for different $$\tau_f$$, when $$T = 40^\circ C$$','interpreter','latex');

%%
% [theta_dot,T] = meshgrid([1e-7:1e-4:1],40:1:80);
% [theta_dot,T] = meshgrid([1e-7:1e-3:0.1 0.1:1e-2:0.9 0.9:1e-1:1],40:1:80);
[theta_dot,T] = meshgrid(1e-5:1e-3:1,40:75);
Fric = zeros(size(theta_dot,1),size(T,2));
for n_x = 1:size(theta_dot,1)
    for n_y = 1:size(T,2)
        Fric(n_x,n_y) = FunFriction(500,theta_dot(n_x,n_y),T(n_x,n_y));
    end
end

figure;
% contour(theta_dot,Fric,T,'showtext','on','LabelSpacing',460);
[b,h] = contour(theta_dot,Fric,T);
clabel(b,h,'labelspacing',380,'fontsize',7);
c = colorbar;
c.Label.String = 'Temperature (¡æ)';
colormap('jet');
xlabel('$$\dot{\theta}(rad/s)$$','interpreter','latex');
ylabel('$$\tau_f(Nm)$$','interpreter','latex');
set(gcf,'color','white');
% xlim([0 1]);
% export_fig Fric_T_curve.eps;
% title('Friction curves for different temperatures $T$, when $$\tau_l = 500Nm$$','interpreter','latex');