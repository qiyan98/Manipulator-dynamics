function tau_f0 = FunFric0(theta_dot)
% Incomplete friction model, only partial parameters are known to the
% controller
% Input: tau_load, motor output torque
%        theta_dot, current angular velocity

sgn_theta_dot = sign(theta_dot);

F_c0 = 196.27; F_ctl = 2.34e-2; F_stl = 0.126; theta_dot_stl = 0.038; alpha = 1.36;
F_s0 = -157.77; F_st = 10.1; theta_dot_s0 = -0.1022; theta_dot_st = 0.004;
F_v0 = 199.14; F_vt = 2022.06; T_v0 = 20.71;
T_0 = 35; tau_load_0 = 1000; % constant
F_c = F_c0 + F_ctl*abs(tau_load_0) + F_stl*abs(tau_load_0)*exp(-abs(theta_dot/theta_dot_stl)^alpha);
F_s = (F_s0 + F_st*T_0)*exp(-abs(theta_dot/(theta_dot_s0 + theta_dot_st*T_0))^alpha);
F_v = (F_v0 + F_vt*exp(-T_0/T_v0))*abs(theta_dot);

% if sgn_theta_dot == 0 % if it is static
%     tau_f0 = min(F_c+F_s+F_v,abs(tau_load))*sign(tau_load);
% elseif  sgn_theta_dot~= 0 % if it is not static
%     tau_f0 = (F_c + F_s + F_v)*sgn_theta_dot;
% end

% tau_f0 = (F_c + F_s + F_v)*sgn_theta_dot;
% tau_f0 = 219.67*sgn_theta_dot + 321*exp(-abs(theta_dot/0.038)^1.36)*sgn_theta_dot + 572.25*theta_dot;
tau_f0 = 219.67*sgn_theta_dot + 572.25*theta_dot;
end