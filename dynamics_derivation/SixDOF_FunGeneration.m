% To generate all necessary anonymous functions used for simulation and
% control in order to simplify codes.

% Trajectory related functions
FunThetaDot = diff(FunTheta,t);
FunThetaDdot = diff(FunThetaDot,t);
FunThetaInt = int(FunTheta);
FunTheta = matlabFunction(FunTheta);
FunThetaDot = matlabFunction(FunThetaDot); %@(t)
FunThetaDdot = matlabFunction(FunThetaDdot); %@(t)
FunThetaInt = matlabFunction(FunThetaInt);% @(t)

% Dynamics related functions
FunM22 = matlabFunction(M22); % @(t3)
FunG2 = matlabFunction(G2); % @(t2,t3)
FunM22_real = matlabFunction(M22_real); % @(t3)
FunG2_real = matlabFunction(G2_real); % @(t2,t3)
FunDesTorque2 = matlabFunction(DesTorque2); %@(t2,t3,t2_ddot)
FunDesTorque2 = @(t,t3)FunDesTorque2(FunTheta(t),t3,FunThetaDdot(t));% @(t,t3)

% State space variables
% X = zeros(2,1); % [x1,x2]' = [t2,t2_dot]' 2x1
% X_dot = zeros(2,1);% [x1_dot,x2_dot]' = [x2,t2_ddot]' 2x1

% DesTorque2_0 = @(t3)FunDesTorque2(0,t3); % scalar
% FunFriction(tau_a,theta_dot,tau_load,T)

% Clear symbolic variables from dynamics
% Note that M22 is constant and time t should be kept.
clear M22 C22 DesTorque2 G2;