function tau_f = FunFriction(tau_load,theta_dot,T)
% Calculate the friciontl torque at given applied torque and velocity
% input: 
% tau_load - output torque, scalar; 
% theta_dot - current joint angular velociy, scalar
% T - temperature, Celsius degree
% output: frictional torque consisting of stick and slip phase

sgn_theta_dot = sign(theta_dot);
% input parameters
% orginial value, velocity is motor side, tau_load is dimensionless
% F_c0 = 3.11e-2; F_ctl = 2.34e-2; F_stl = 0.126; theta_dot_stl = 9.22; alpha = 1.36;
% F_s0 = -2.5e-2; F_st = 1.6e-3; theta_dot_s0 = -24.81; theta_dot_st = 0.98;
% F_v0 = 1.3e-4; F_vt = 1.32e-3; T_v0 = 20.71; rv_i = 242.73; max_torque = 26;

% modified value, velocity is joint side, tau_load is actual load
% Max torque = 6310.98 Nm, rv_i = 242.73;
F_c0 = 196.27; F_ctl = 2.34e-2; F_stl = 0.126; theta_dot_stl = 0.038; alpha = 1.36;
F_s0 = -157.77; F_st = 10.1; theta_dot_s0 = -0.1022; theta_dot_st = 0.004;
F_v0 = 199.14; F_vt = 2022.06; T_v0 = 20.71;
F_c = F_c0 + F_ctl*abs(tau_load) + F_stl*abs(tau_load)*exp(-abs(theta_dot/theta_dot_stl)^alpha);
F_s = (F_s0 + F_st*T)*exp(-abs(theta_dot/(theta_dot_s0 + theta_dot_st*T))^alpha);
F_v = (F_v0 + F_vt*exp(-T/T_v0))*abs(theta_dot);
% The above F_c, F_s, F_v is absolute value, irrespective of direction.

if sgn_theta_dot == 0 % if it is static
    tau_f = min(F_c+F_s+F_v,abs(tau_load))*sign(tau_load);
elseif  sgn_theta_dot~= 0 % if it is not static
    tau_f = (F_c + F_s + F_v)*sgn_theta_dot;
end
end
