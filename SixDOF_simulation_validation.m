function SixDOF_simulation_validation(FunTheta,MaxTime)
% Validation for simulation, comparison between different relative
% tolerance is made.
% Input: FunTheta and MaxTime
Tol = [1e-3 1e-5 1e-7 1e-8];
cur_t3 = 0; T = 40;

SixDOF_dynamics_simplified; % load necessary dynamic variables
SixDOF_FunGeneration; % load necessary anonymous functions

FunDesTorque2 = @(t)FunDesTorque2(t,cur_t3);
M22 = FunM22(cur_t3);
FunG2 = @(t2) FunG2(t2,cur_t3);

% kp = 200*1e2; kd = 50*1e2; ki = 5*1e2;
% error_pos = @(t,X)( FunTheta(t) - X(1));
% error_vel = @(t,X)( FunThetaDot(t) - X(2));
% error_int = @(t,X)( FunThetaInt(t) - X(3));
% FunOutput = @(t,X) 1*(kp*error_pos(t,X) + kd*error_vel(t,X) + ki*error_int(t,X));
% X_0 = [FunTheta(0);FunThetaDot(0);FunThetaInt(0)];
% SixDOF_Simplified_EOM = @(t,X)[X(2);M22^-1*FunU(t,X);X(1)];

FunOutput = @(t,X) FunDesTorque2(t);
X_0 = [FunTheta(0);FunThetaDot(0)];
FunU = @(t,X) FunOutput(t,X) - FunG2(X(1)) - FunFriction(FunOutput(t,X),X(2),T);
SixDOF_Simplified_EOM = @(t,X)[X(2);M22^-1*FunU(t,X)];

% solve dynamic equation and obtain necessary data for visualization
data_cell = cell(length(Tol),3); str_cell = cell(1,length(Tol));
for j = 1:length(Tol)
    opts = odeset('RelTol',Tol(j),'refine',10);
    [t,X] = ode45(SixDOF_Simplified_EOM,[0 MaxTime],X_0,opts);
    data_cell{j,1} = t; data_cell{j,2} = X; data_cell{j,3} = gradient(X(:,2),t);
    clear t X;
    str_cell{1,j} = ['Tolerance = ' num2str(Tol(j))];
end

%%
figure; hold on;
cell_linespec = cell(6,1);
cell_linespec{1} = '-o'; cell_linespec{2} = '-*';
cell_linespec{3} = '-d'; cell_linespec{4} = '-^'; 
cell_linespec{5} = '-.>'; cell_linespec{6} = '-.h';
for j = 1:length(Tol)
    plot(data_cell{j,1},data_cell{j,3},cell_linespec{j},'MarkerIndices',1:4:length(data_cell{j,3}));
%     plot(data_cell{j,1},data_cell{j,2}(:,2));
end
xlim([2.15,2.55]);
box on;
xlabel('Time (s)');
ylabel('Acceleration (rad/s^2)');
legend(str_cell);
set(gcf, 'Color', 'w');
export_fig computation_validation.eps;
% tightfigadv;
end