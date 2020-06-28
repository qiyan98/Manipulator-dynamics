% Present dynamic equation from computed results
% This file presents simplified 1-DOF (for joint #2) situation for ABB IRB 6620.
syms t t2 t2_dot t2_ddot t3 real;
syms a2 real;
syms g m2 m3 xc2 yc2 zc2 xc3 yc3 zc3 real;
syms Ixx2 Iyy2 Izz2 Ixy2 Iyz2 Ixz2 real;
syms Ixx3 Iyy3 Izz3 Ixy3 Iyz3 Ixz3 real;

% M22 = m2*xc2^2 + m2*yc2^2 + Izz2;
% C22 = 0;
% G2 = -g*m2*(xc2*cos(t2) - yc2*sin(t2));

M22 = m3*a2^2 + 2*m3*cos(t3)*a2*xc3 - 2*m3*sin(t3)*a2*yc3 + m2*xc2^2 + m3*xc3^2 + m2*yc2^2 + m3*yc3^2 + Izz2 + Izz3;
G2 = g*(a2*m3*cos(t2) + m2*xc2*cos(t2) - m2*yc2*sin(t2) + m3*xc3*cos(t2 - t3) + m3*yc3*sin(t2 - t3));
M22_real = M22;
G2_real = G2;

M22 = subs(M22,[m2 xc2 yc2 zc2],[110.698 0.418 0 -0.218]);
M22 = subs(M22,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[0.011 -0.005 -11.122 6.208 35.317 29.285]);
M22 = subs(M22,[m3 xc3 yc3 zc3],[250.044 0.156 0.268 -0.001]);
M22 = subs(M22,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456]);
M22 = subs(M22,a2,0.975);
M22 = eval(M22);

p_para = 1.1;
p_inertia_para = 0.9;
M22_real = subs(M22_real,[m2 xc2 yc2 zc2],[p_para*110.698 0.418 0 -0.218]);
M22_real = subs(M22_real,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[0.011 -0.005 -11.122 6.208 35.317 29.285*p_inertia_para]);
M22_real = subs(M22_real,[m3 xc3 yc3 zc3],[p_para*250.044 0.156 0.268 -0.001]);
M22_real = subs(M22_real,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456*p_inertia_para]);
M22_real = subs(M22_real,a2,0.975);
M22_real = eval(M22_real);

G2 = subs(G2,[m2 xc2 yc2 zc2],[110.698 0.418 0 -0.218]);
G2 = subs(G2,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[0.011 -0.005 -11.122 6.208 35.317 29.285]);
G2 = subs(G2,[m3 xc3 yc3 zc3],[250.044 0.156 0.268 -0.001]);
G2 = subs(G2,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456]);
G2 = subs(G2,[a2 g],[0.975 9.81]);

G2_real = subs(G2_real,[m2 xc2 yc2 zc2],[p_para*110.698 0.418 0 -0.218]);
G2_real = subs(G2_real,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[0.011 -0.005 -11.122 6.208 35.317 29.285*p_inertia_para]);
G2_real = subs(G2_real,[m3 xc3 yc3 zc3],[p_para*250.044 0.156 0.268 -0.001]);
G2_real = subs(G2_real,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456*p_inertia_para]);
G2_real = subs(G2_real,[a2 g],[0.975 9.81]);

DesTorque2 = M22*t2_ddot + G2;

clear t2 t2_dot t2_ddot t3 a2;
clear g m2 m3 xc2 yc2 zc2 xc3 yc3 zc3;
clear Ixx2 Iyy2 Izz2 Ixy2 Iyz2 Ixz2;
clear Ixx3 Iyy3 Izz3 Ixy3 Iyz3 Ixz3;