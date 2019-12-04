clear;
close all;
% Numerical validation for ode solution for manipulator dynamics simulation
% A 1-DOF sliding mass-spring system with friction is adopted from Marques,
% 2016

% X = [x1 x2] = [x x_dot]
% X' = [x1' x2'] = [x2 x_ddot]
k = 2; % spring stiffness
m = 1; % mass of the block
v_b = 0.1; % velocity of belt
FunU = @(t,X) -Validation_fun_friction(X(2) - v_b,-k*X(1)) -k*X(1);
EOM = @(t,X) [X(2);m^-1*FunU(t,X)]; % equation of motion
X_0 = [0 0.1];

% ODE solution process
MaxTime = 20;
opts = odeset('RelTol',1e-6);
[t,X] = ode45(EOM,[0 MaxTime],X_0,opts);
% opts = odeset('RelTol',1e-6);
% [t,X] = ode45(EOM,[0:1e-3:MaxTime],X_0);

%% Visulization of results
% figure;
% v = [-3:0.1:-0.1,-0.1:1e-5:0.1,0.1:0.1:3]; fric = zeros(10,1);
% for j = 1:length(v)
%     for f_a = -2:0.01:2
%         fric(j) = Validation_fun_friction(v(j),f_a);
%     end
% end
% plot(v,fric);
% xlabel('Velocity(m/s)');
% ylabel('Friction force(N)');

width = 1.5;
figure;
plot(t,X(:,1),'linewidth',width);
ylim([0 1.5]);
xlabel('time(s)');
ylabel('displacement(m)');

figure;
plot(t,-(X(:,2)-v_b),'linewidth',width);
ylim([0 1]);
xlabel('time(s)');
ylabel('Relative velocity (m/s)');

figure;
fric = zeros(10,1);
for j = 1:length(t)
    fric(j) = -Validation_fun_friction(X(j,2) - v_b,-k*X(j,1));
end
plot(t,fric,'linewidth',width);
ylim([0 2.5]);
% hold on;
% plot(t,-k*X(:,1));
xlabel('time(s)');
ylabel('Friction force(N)');


% figure;
% x_ddot = gradient(X(:,2));
% plot(t,x_ddot);
% xlabel('time(s)');
% ylabel('Relative velocity (m/s^2)');