function force_f = Validation_fun_friction(v,f_a)
% Input: velocity of the mass block with regard to ground
%        external force exerted on the mass block
% Output: friction force between conveyor belt and mass block

% Set motion-related parameters
m = 1; g = 9.81;

% Set frictional parameters from Marques, 2016
f_c = 0.1*m*g; f_s = 0.15*m*g; v_s = 0.001; f_v = 0.1;

if v == 0
    force_f = min(abs(f_a),f_s)*sign(f_a);
else
    force_f = (f_c + (f_s - f_c)*exp(-(abs(v)/v_s)^2))*sign(v);
    force_f = force_f + f_v*v;
end
end