% Take ABB IRB 6620 as example for simulation by manual Newton-Euler
% Recursive methods
% Parameters are taken from Sergi and Daniel, 2011 with standard D-H chain.
close all;
clear;
SixDOF_symbolics;
%% Kinematics - frame translation
T10 = [cos(t1) -sin(t1) 0 0;
    sin(t1) cos(t1) 0 0;
    0 0 1 d1;
    0 0 0 1];
T21 = [cos(t2) -sin(t2) 0 a1;
    0 0 -1 0;
    sin(t2) cos(t2) 0 0;        
    0 0 0 1];
T32 = [cos(t3) -sin(t3) 0 a2;
       -sin(t3) -cos(t3) 0 0;
       0 0 -1 0;
       0 0 0 1];
T30 = T10*T21*T32;
T30 = simplify(T30);
R10 = T10(1:3,1:3); R01 = R10^-1;   
R21 = T21(1:3,1:3); R12 = R21^-1;
R32 = T32(1:3,1:3); R23 = R32^-1;
Fun1 = matlabFunction(T30);
%% Kinematics - velocity & acceleration propagation

Omega = cell(1,3); % cell for omega of joints
OmegaAcc = cell(1,3); % cell for angular acceleration
Vel = cell(1,3); % cell for linear velocity of joints
VelCen = cell(1,3); % cell for linear volocity of mass centers of links
VelAcc = cell(1,3); % cell for joint acceleration
VelCenAcc = cell(1,3); % cell for COG acceleration

Vel_0 = [0 0 0]'; Omega_0 = [0 0 0]';
OmegaAcc_0 = [0 0 0]'; VelAcc_0 = [0 0 g]';
% sequence: position, next omega, next joint velocity, next COG velociy ...
% next angular acceleration, next joint acceleration
P10 = [0 0 d1]';
Pc1 = [xc1 yc1 zc1]';
Omega{1} = R01*Omega_0 + [0 0 t1_dot]';
Vel{1} = R01*(Vel_0 + cross(Omega_0,P10));
VelCen{1} = Vel{1} + cross(Omega{1},Pc1);
OmegaAcc{1} = R01*OmegaAcc_0 + cross(R01*Omega_0,[0 0 t1_dot]') + [0 0 t1_ddot]';
VelAcc{1} = R01*(cross(OmegaAcc_0,P10) + cross(Omega_0,cross(Omega_0,P10)) + VelAcc_0);
VelCenAcc{1} = VelAcc{1} + cross(Omega{1},cross(Omega{1},Pc1)) + cross(OmegaAcc{1},Pc1);

P21 = [a1 0 0]';
Pc2 = [xc2 yc2 zc2]';
Omega{2} = R12*Omega{1} + [0 0 t2_dot]';
Vel{2} = R12*(Vel{1} + cross(Omega{1},P21));
VelCen{2} = Vel{2} + cross(Omega{2},Pc2);
OmegaAcc{2} = R12*OmegaAcc{1} + cross(R12*Omega{1},[0 0 t2_dot]') + [0 0 t2_ddot]';
VelAcc{2} = R12*(cross(OmegaAcc{1},P21) + cross(Omega{1},cross(Omega{1},P21)) + VelAcc{1});
VelCenAcc{2} = VelAcc{2} + cross(Omega{2},cross(Omega{2},Pc2)) + cross(OmegaAcc{2},Pc2);

P32 = [a2 0 0]';
Pc3 = [xc3 yc3 zc3]';
Omega{3} = R23*Omega{2} + [0 0 t3_dot]';
Vel{3} = R23*(Vel{2} + cross(Omega{2},P32));
VelCen{3} = Vel{3} + cross(Omega{3},Pc3);
OmegaAcc{3} = R23*OmegaAcc{2} + cross(R23*Omega{2},[0 0 t3_dot]') + [0 0 t3_ddot]';
VelAcc{3} = R23*(cross(OmegaAcc{2},P32) + cross(Omega{2},cross(Omega{2},P32)) + VelAcc{2});
VelCenAcc{3} = VelAcc{3} + cross(Omega{3},cross(Omega{3},Pc3)) + cross(OmegaAcc{3},Pc3);
%% Symbolic Dynamcis using Newton-Euler Recursive method
% Outward iterations
Force = cell(1,3);
Moment = cell(1,3);
Force{1} = m1*VelCenAcc{1}; Force{2} = m2*VelCenAcc{2}; Force{3} = m3*VelCenAcc{3};
Moment{1} = I1*OmegaAcc{1} + cross(Omega{1},I1*Omega{1});
Moment{2} = I2*OmegaAcc{2} + cross(Omega{2},I2*Omega{2});
Moment{3} = I3*OmegaAcc{3} + cross(Omega{3},I3*Omega{3});
% Inward iterations
force = cell(1,3); moment = cell(1,3);
force{3} = Force{3};
moment{3} = Moment{3} + cross(Pc3,Force{3});
force{2} = Force{2} + R32*force{3};
moment{2} = Moment{2} + R32*moment{3} + cross(Pc2,Force{2}) + cross(P32,R32*force{3});
force{1} = Force{1} + R21*force{2};
moment{1} = Moment{1} + R21*moment{2} + cross(Pc1,Force{1}) + cross(P21,R21*force{2});
% Evaluated torque
torque = simplify([moment{1}(3) moment{2}(3) moment{3}(3)]');
torque_23 = subs(torque,[t1 t1_dot t1_ddot t3_dot t3_ddot],[0 0 0 0 0]);

Torq = torque;
% Torq = subs(Torq,[m1 xc1 yc1 zc1],[220.489 0.419 0.033 -0.151]);
% Torq = subs(Torq,[Ixy1 Iyz1 Ixz1 Ixx1 Iyy1 Izz1],[1.826 -0.176 -0.543 14.744 25.283 17.320]);
Torq = subs(Torq,[m2 xc2 yc2 zc2],[110.698 0.418 0 0.218]);
Torq = subs(Torq,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[-0.011 -0.005 11.122 6.208 35.317 30.285]);
Torq = subs(Torq,[m3 xc3 yc3 zc3],[250.440 0.156 0.268 -0.001]);
Torq = subs(Torq,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456]);
Torq = subs(Torq,[a1 d1 a2 g],[0.32 0.68 0.975 9.81]);
Torq = subs(Torq,[t1 t2 t3],[0 50 50]./57.3);
Torq = subs(Torq,[t1_dot t2_dot t3_dot],[10 10 10]./57.3);
Torq = subs(Torq,[t1_ddot t2_ddot t3_ddot],[10 10 10]./57.3);
vpa(Torq,5)
Fun1 = @(t1,t2,t3)Fun1(0.32,0.975,0.68,t1,t2,t3);

% Input these manually
M22 = m3*a2^2 + 2*m3*cos(t3)*a2*xc3 - 2*m3*sin(t3)*a2*yc3 + m2*xc2^2 + m3*xc3^2 + m2*yc2^2 + m3*yc3^2 + Izz2 + Izz3;
G2 = g*(a2*m3*cos(t2) + m2*xc2*cos(t2) - m2*yc2*sin(t2) + m3*xc3*cos(t2 - t3) + m3*yc3*sin(t2 - t3));
tau_2 = M22*t2_ddot + G2;
error_2 = simplify(torque_23(2) - tau_2);

% M11 = Ixx2 + Izz1 + a1^2*m2 + m1*xc1^2 + m1*yc1^2 + m2*yc2^2 + m2*zc2^2 - Ixx2*cos(t2)^2 + Iyy2*cos(t2)^2 - Ixy2*sin(2*t2) + m2*xc2^2*cos(t2)^2 - m2*yc2^2*cos(t2)^2 - m2*xc2*yc2*sin(2*t2) + 2*a1*m2*xc2*cos(t2) - 2*a1*m2*yc2*sin(t2);
% M12 = Iyz2*cos(t2) + Ixz2*sin(t2) + m2*yc2*zc2*cos(t2) + m2*xc2*zc2*sin(t2);
% M21 = Iyz2*cos(t2) + Ixz2*sin(t2) + m2*yc2*zc2*cos(t2) + m2*xc2*zc2*sin(t2);
% M22 = m2*xc2^2 + m2*yc2^2 + Izz2;