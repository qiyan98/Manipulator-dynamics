function SixDOF_simulation(CurT,Configuration,FunTheta,MaxTime,mode_ctrl,options)
% Numerical Dynamic simulation - dynamic simulation with complex friction
% ODEs solved by ode45 or customized discrete solver, for joint #2 in ABB IRB 6620 only.
% Open loop simulation is done and its resutls is present, edited based on
% SixDOF_TrajSim2.m

% Input: parameters that could be specified
% CurT - a vector of current temperature
% Configuration - a vector of t3 angles
% FunTheta - an anonymous function for set trajectory
% MaxTime - maximum simulation time
% options - options to set how resutls is presented


% Deal with input parameters/arguments
if nargin < 6
    options = 'true';
end

if ~strcmp(options,'true') && ~strcmp(options,'error')
    disp('Wrong command for options! True values would be presented.');
    return;
    options = 'true';
end

% if ~strcmp(mode_ctrl,'openloop') && ~strcmp(mode_ctrl,'ctc') && ~strcmp(mode_ctrl,'pid')
%     disp('Wrong command for mode! Open loop simulation would be conducted.');
%     mode_ctrl = 'openloop';
% end

if strcmpi(mode_ctrl,'openloop')
    mode_ctrl = 0;
elseif strcmpi(mode_ctrl,'pid')
    mode_ctrl = 1;
elseif strcmpi(mode_ctrl,'ctc')
    mode_ctrl = 2;
elseif strcmpi(mode_ctrl,'FNTSM')
    mode_ctrl = 3;
elseif strcmpi(mode_ctrl,'AFNTSM_SF')
    mode_ctrl = 4;
elseif strcmpi(mode_ctrl,'AFNTSM_LT')
    mode_ctrl = 5;
elseif strcmpi(mode_ctrl,'LSMC')
    mode_ctrl = 6;
elseif strcmpi(mode_ctrl,'hinf')
    mode_ctrl = 7;
else
    mode_ctrl = -1;
    disp([mode_ctrl,' is wrong command for control mode!']);
    save SixDOF_simData.mat mode_ctrl;
    return;
end

MaxTime = max(MaxTime,1); % The time should not be less than 1 second.

% Prepare for ode solver
SixDOF_dynamics_simplified; % load necessary dynamic variables
SixDOF_FunGeneration; % load necessary anonymous functions
% create variables to store calculation results
n_t3 = length(Configuration);
n_tem = length(CurT);
% values dependent on mode
if mode_ctrl == 0 % open loop condition
    X_0 = [FunTheta(0);FunThetaDot(0)];
    data = cell(n_t3,n_tem,4);
elseif mode_ctrl == 1 % pid condition
    X_0 = [FunTheta(0);FunThetaDot(0);FunThetaInt(0)];
    data = cell(n_t3,n_tem,5);
elseif mode_ctrl == 2 % ctc condition
    X_0 = [FunTheta(0);FunThetaDot(0)];
    data = cell(n_t3,n_tem,5);
elseif mode_ctrl == 3 % FNTSM condition
    X_0 = [FunTheta(0);FunThetaDot(0)];
    data = cell(n_t3,n_tem,5);
elseif mode_ctrl == 4 % AFNTSM_SF condition
    X_0 = [FunTheta(0);FunThetaDot(0);0;0;0];
    data = cell(n_t3,n_tem,5);
elseif mode_ctrl == 5 % AFNTSM-TL condition
    X_0 = [FunTheta(0);FunThetaDot(0);0;0;0];
    data = cell(n_t3,n_tem,7);
elseif mode_ctrl == 6 % LSMC condition
    X_0 = [FunTheta(0);FunThetaDot(0)];
    data = cell(n_t3,n_tem,5);
elseif mode_ctrl == 7 % H-infinity condition
    X_0 = [FunTheta(0);FunThetaDot(0)];
    data = cell(n_t3,n_tem,5);
end
% data illustration:
% 1st dimension - different configuration
% 2nd dimension - different temperature
% 3rd dimension - t-X-t2_ddot-FricTorq

S = RandStream('mt19937ar','Seed',5489); % randstream for awgn

% Solve the equation of motion by ode45 and collect its results
FunDesTorque2_backup = FunDesTorque2;
FunM22_backup = FunM22;
FunG2_backup = FunG2;
FunM22_real_backup = FunM22_real;
FunG2_real_backup = FunG2_real;
for m = 1:n_t3 % index of configurations
    FunDesTorque2 = @(t)FunDesTorque2_backup(t,Configuration(m));
    M22 = FunM22_backup(Configuration(m));
    FunG2 = @(t2) FunG2_backup(t2,Configuration(m));
    for j = 1:n_tem % index of temperatures
        clear t X;
        p_d_bar = 200; % disturbance parameters
        dt_t3 = 1; vel_t3 = -70/57.3;
        FunVary_t3 = @(t) Configuration(m) + (t>4 & t<4+dt_t3)*(t-4)*vel_t3 ...
            + (t>=4+1)*dt_t3*vel_t3;
        FunVaryingT = @(t) CurT(j) + (t>6.8)*20;
        
        FunVaryingM22_real = @(t) FunM22_real_backup(FunVary_t3(t));
        FunVaryingG2_real = @(t,X) FunG2_real_backup(X(1),FunVary_t3(t));
        FunVaryingM220 = @(t) FunM22_backup(FunVary_t3(t));
        FunVaryingG20 = @(t,X) FunG2_backup(X(1),FunVary_t3(t));
        
        %%%%% shared components for all sliding mode controls %%%%%
        p_lambda = 0.01; p_gamma = 1.3; p_rho = 0.9; p_delta = 0.05;
        Fun_k1 = @(t,X) 80; Fun_k2 = @(t,X) 70;
        mu_0 = 20; mu_1 = 20; mu_2 = 20;
        
        snr_theta = 20; % fixed
        snr_thetaDot = 20; % fixed
        t_noise = 0.0;
        
        error_pos = @(t,X)(awgn(X(1),snr_theta,'measured',S) - FunTheta(t)+...
            0*FunNoise(t,1)); % q - q_r, agwn
        error_vel = @(t,X)(awgn(X(2),snr_thetaDot,'measured',S) - FunThetaDot(t)+...
            0*FunNoise(t,2)); % q_dot - q_dot_r, awgn
        
        Fun_s = @(t,X) error_pos(t,X) + p_lambda*signum(error_vel(t,X),p_gamma);
        %%%%% shared components for all sliding mode controls %%%%%
        
        %%%%% components for AFNTSM %%%%%
        % a_0_hat - X(3); a_1_hat - X(4); a_2_hat - X(5)
        Fun_a0d = @(t,X) mu_0*abs(error_vel(t,X))^(p_gamma-1)*abs(Fun_s(t,X));
        Fun_a1d = @(t,X) mu_1*abs(error_vel(t,X))^(p_gamma-1)*...
            abs(awgn(X(1),snr_theta,'measured',S))*abs(Fun_s(t,X));
        Fun_a2d = @(t,X) mu_2*abs(error_vel(t,X))^(p_gamma-1)*...
            abs(awgn(X(2),snr_thetaDot,'measured',S))^2*abs(Fun_s(t,X));
        Fun_xi = @(t,X) X(3)+X(4)*abs(awgn(X(1),snr_theta,'measured',S))...
            + X(5)*abs(awgn(X(2),snr_thetaDot,'measured',S))^2; % adaptive coefficient
        %%%%% components for AFNTSM %%%%%
        
        FunD = @(t) p_d_bar*sin(t) + p_d_bar/10*sin(200*pi*t) + (t>7.4&t<(7.4+0.1))*10*p_d_bar; % uncertainty disturbance

        if mode_ctrl == 0 % open loop
            FunOutput = @(t,X) FunDesTorque2(t);
        elseif mode_ctrl == 1 % pid
            kp = 5e5; kd = 8.5e4; ki = 2e4;
            error_pos = @(t,X)( FunTheta(t) - X(1));
            error_vel = @(t,X)( FunThetaDot(t) - X(2));
            error_int = @(t,X)( FunThetaInt(t) - X(3));
            FunOutput = @(t,X) 1*(kp*error_pos(t,X) + kd*error_vel(t,X) + ki*error_int(t,X));
        elseif mode_ctrl == 2 % ctc
            kp = 300; kd = 70;
            error_pos = @(t,X) -error_pos(t,X);
            error_vel = @(t,X) -error_vel(t,X);
            FunTau_pd = @(t,X) FunVaryingM220(t)*(kp*error_pos(t,X) + kd*error_vel(t,X));
            Fun_dyn = @(t,X) FunVaryingM220(t)*FunThetaDdot(t) + FunTau_pd(t,X) + FunVaryingG20(t,X);
            FunEquilibrium = @(t,X,y) Fun_dyn(t,X) + FunFriction(y,X(2),FunVaryingT(t)) - y;
            FunOutput = @(t,X)fzero(@(y)FunEquilibrium(t,X,y),FunVaryingM220(t)*FunThetaDdot(t) + FunVaryingG20(t,X)); % solve for the implicit equation for output torque!
        elseif mode_ctrl == 3 % FNTSM
            FunTau0_backup = @(t,X) FunVaryingM220(t)*FunThetaDdot(t) + FunVaryingG20(t,X) - FunVaryingM220(t)/(p_lambda*p_gamma)*signum(error_vel(t,X),2-p_gamma);
            FunTau1 = @(t,X) -FunVaryingM220(t)*(Fun_k1(t,X)*Fun_s(t,X) + Fun_k2(t,X)*signum(Fun_s(t,X),p_rho));
            FunOutput = @(t,X) FunTau0_backup(t,X) + FunFric0(X(2)) + FunTau1(t,X);
        elseif mode_ctrl == 4 % AFNTSM-SF
            FunTau0_backup = @(t,X) FunVaryingM220(t)*FunThetaDdot(t) + FunVaryingG20(t,X) - FunVaryingM220(t)/(p_lambda*p_gamma)*signum(error_vel(t,X),2-p_gamma);
            FunTau1 = @(t,X) -FunVaryingM220(t)*(Fun_k1(t,X)*Fun_s(t,X) + Fun_k2(t,X)*signum(Fun_s(t,X),p_rho) + Fun_xi(t,X)*FunSat(Fun_s(t,X)/p_delta));
            FunOutput = @(t,X) FunTau0_backup(t,X) + FunFric0(X(2)) + FunTau1(t,X);
        elseif mode_ctrl == 5 % AFNTSM-TL
            FunTau0_backup = @(t,X) FunVaryingM220(t)*FunThetaDdot(t) + FunVaryingG20(t,X) - FunVaryingM220(t)/(p_lambda*p_gamma)*signum(error_vel(t,X),2-p_gamma);
            FunTau1 = @(t,X) -FunVaryingM220(t)*(Fun_k1(t,X)*Fun_s(t,X) + Fun_k2(t,X)*signum(Fun_s(t,X),p_rho) + Fun_xi(t,X)*FunSat(Fun_s(t,X)/p_delta));
            FunEquilibrium = @(t,X,y) FunTau0_backup(t,X) + FunFriction(y,X(2),FunVaryingT(t)) - y + FunTau1(t,X);
            FunOutput = @(t,X)fzero(@(y)FunEquilibrium(t,X,y),FunTau0_backup(t,X) + FunFric0(X(2)) + FunTau1(t,X));
        elseif mode_ctrl == 6 % LSMC
            p_c = 5; p_k = 2000;
            Fun_s = @(t,X) p_c*error_pos(t,X) + error_vel(t,X);
            FunTau0_backup = @(t,X) FunVaryingM220(t)*(FunThetaDdot(t) - p_c*error_vel(t,X)) + FunVaryingG20(t,X);
            FunTau1 = @(t,X) -p_k*FunSat(Fun_s(t,X)/p_delta);
            FunEquilibrium = @(t,X,y) FunTau0_backup(t,X) + FunFriction(y,X(2),FunVaryingT(t)) - y + FunTau1(t,X);
            FunOutput = @(t,X)fzero(@(y)FunEquilibrium(t,X,y),FunTau0_backup(t,X) + FunFric0(X(2)) + FunTau1(t,X));
        elseif mode_ctrl == 7 % H-infinity control
            p_alpha1 = 1; p_alpha2 = 1; p_rho = 0.005; p_xi = 0.1;
            A_z = [0 1;-p_alpha2 -p_alpha1]; B_z = [0;1];
            Q_z = 5*eye(2) + p_xi*eye(2);
            [P,~,~] = care(A_z,B_z/p_rho,Q_z);
            K = 1/p_rho^2*B_z'*P;
            FunTau0 = @(t,X) FunVaryingG20(t,X) + 572.25*(X(2));
            Fun_z = @(t,X) [error_pos(t,X); error_vel(t,X)];
            FunTau1 = @(t,X) FunVaryingM220(t)*(FunThetaDdot(t) -p_alpha1*error_vel(t,X) -p_alpha2*error_pos(t,X) -K*Fun_z(t,X) );
            FunOutput = @(t,X) FunTau0(t,X) + FunTau1(t,X);
        end
        
        FunOutput = @(t,X) 6311*FunSat(FunOutput(t,X)/6311);
        FunU = @(t,X) FunOutput(t,X) - FunVaryingG2_real(t,X) - FunFriction(FunOutput(t,X),X(2),FunVaryingT(t));
        FunU = @(t,X) FunU(t,X) - FunD(t);
        if mode_ctrl == 1 % pid, 3rd order system
            SixDOF_Simplified_EOM = @(t,X)[X(2);FunVaryingM22_real(t)^-1*FunU(t,X);X(1)];
        elseif mode_ctrl == 4||mode_ctrl == 5 % AFNTSM-SF/AFNTSM-LT
            SixDOF_Simplified_EOM = @(t,X)[X(2);FunVaryingM22_real(t)^-1*FunU(t,X);Fun_a0d(t,X);Fun_a1d(t,X);Fun_a2d(t,X)];
        else % not pid, 2nd system :openloop/ctc/FNTSM/H-inf
            SixDOF_Simplified_EOM = @(t,X)[X(2);FunVaryingM22_real(t)^-1*FunU(t,X)];
        end

        opts = odeset('RelTol', 1e-3,'AbsTol',1e-6,'Refine',6); % debug
%         [t_ode,X_ode] = ode45(SixDOF_Simplified_EOM,[0,MaxTime],X_0,opts); % debug

        t_ode = 0:1e-3:MaxTime;
        opts = sdeset('OutputFUN',@sdeplot,...
              'SDEType','Stratonovich',...
              'RandSeed',2);
        if mode_ctrl == 1 % pid, 3rd order system
            X_ode = sde_euler(SixDOF_Simplified_EOM,@(t,y)[0,0,0]',t_ode...
            ,X_0,opts);
        elseif mode_ctrl == 4||mode_ctrl == 5 % AFNTSM-SF/AFNTSM-LT
            X_0t = [FunTheta(t_ode(1));FunThetaDot(t_ode(1));0;0;0];
            X_ode = sde_euler(SixDOF_Simplified_EOM,@(t,y)[0,0,0,0,0]',t_ode...
            ,X_0t,opts);
        else % not pid, 2nd system :openloop/ctc/FNTSM/H-inf
            X_ode = sde_euler(SixDOF_Simplified_EOM,@(t,y)[0,0]',t_ode...
            ,X_0,opts);
        end
        
        X_ode = sde_milstein(SixDOF_Simplified_EOM,...
                [0,0,0,0,0]',t_ode,X_0,opts);
            
        T2Ddot = gradient(X_ode(:,2),t_ode);
        FricTorq2 = zeros(10,1); % actual friction torque
        FricDist = zeros(10,1);  % frictional disturbance
        FircDistT = zeros(10,1); % frictional disturbance caused by temperature
        FricDistL = zeros(10,1); % frictional disturbance caused by load
        RealTorq2 = zeros(10,1); % actual output torque
        for n = 1:length(t_ode)
            RealTorq2(n) = FunOutput(t_ode(n),X_ode(n,:));
            FricTorq2(n) = FunFriction(RealTorq2(n),X_ode(n,2),FunVaryingT(t_ode(n)));
            FricDist(n) = FricTorq2(n) - FunFriction(1000,X_ode(n,2),35); % 35 C,1k Nm
            FircDistT(n) = FricTorq2(n) - FunFriction(RealTorq2(n),X_ode(n,2),35);
            FricDistL(n) = FricTorq2(n) - FunFriction(1000,X_ode(n,2),FunVaryingT(t_ode(n)));
        end
        data{m,j,1} = t_ode; data{m,j,2} = X_ode;
        data{m,j,3} = T2Ddot; data{m,j,4} = FricTorq2;
        if mode_ctrl ~= 0 % not open loop, i.e. control mode
            data{m,j,5} = RealTorq2;
            data{m,j,6} = FricDist; % frictional disturbance;
            data{m,j,7} = FircDistT; % frictional disturbance;
            data{m,j,8} = FricDistL; % frictional disturbance;
        end
        clear t X T2Dot FricTorq2 RealTorq2;
    end % end of this temperature
end % end of this theta_3

%% Calculate data for error comparison
control_error = cell(n_t3,n_tem,3);
true_traj = cell(n_t3,n_tem,3);
% 1st dimension - different configuration
% 2nd dimension - different temperature
% 3rd dimension - error of position-velociy-acceleration respectively
for m = 1:n_t3
    for j = 1:n_tem
        ThetaVec = zeros(10,1); ThetaDotVec = zeros(10,1); ThetaDdotVec = zeros(10,1);
        for n = 1:length(data{m,j,1})
            ThetaVec(n) = FunTheta(data{m,j,1}(n));
            ThetaDotVec(n) = FunThetaDot(data{m,j,1}(n));
            ThetaDdotVec(n) = FunThetaDdot(data{m,j,1}(n));
        end
        control_error{m,j,1} = ThetaVec - data{m,j,2}(:,1);     % position error

        control_error{m,j,2} = ThetaDotVec - data{m,j,2}(:,2);  % velocity error
        control_error{m,j,3} = ThetaDdotVec - data{m,j,3};      % acceleration error
        true_traj{m,j,1} = ThetaVec;
        true_traj{m,j,2} = ThetaDotVec;
        true_traj{m,j,3} = ThetaDdotVec;
        clear ThetaVec ThetaDotVec ThetaDdotVec;
    end
end
if mode_ctrl == 1 % pid
    save SixDOF_simData_pid.mat;
elseif mode_ctrl == 2 %ctc
    save SixDOF_simData_ctc.mat;
elseif mode_ctrl == 3 %fntsm
    save SixDOF_simData_fntsm.mat;
elseif mode_ctrl == 4 %afntsm
    save SixDOF_simData_afntsm_sf.mat;
elseif mode_ctrl == 5 %afntsm_lt
    save SixDOF_simData_afntsm_lt.mat
elseif mode_ctrl == 6 %lsmc
    save SixDOF_simData_lsmc.mat
elseif mode_ctrl == 7 %h-infinity
    save SixDOF_simData_hinf.mat
end
save SixDOF_simData.mat;

end