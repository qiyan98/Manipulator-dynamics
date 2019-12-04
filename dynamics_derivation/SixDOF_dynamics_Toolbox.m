% Take ABB IRB 6620 as example for simulation by MATLAB Robotics Toolbox
% Parameters are taken from Sergi and Daniel, 2011 with standard D-H chain.
clear;
close all;
SixDOF_symbolics; % symbolic variables for dynamics model

% Version from Sergi and Daniel, 2011
% R1 = Link('revolute','d',680,'a',320,'alpha',pi/2,'m',m1,'r',[xc1 yc1 0]','I',I1);
% R2 = Link('revolute','d',0,'a',975,'alpha',0,'m',m2,'r',[xc2 yc2 0]','I',I2);
% R3 = Link('revolute','d',0,'a',200,'alpha',pi/2,'m',m3,'r',[xc3 yc3 0]','I',I3);
% R4 = Link('revolute','d',887,'a',0,'alpha',-pi/2,'m',m4,'r',[xc4 yc4 0]','I',I4);
% R5 = Link('revolute','d',0,'a',0,'alpha',pi/2,'m',m5,'r',[xc5 yc5 0]','I',I5);
% R6 = Link('revolute','d',300,'a',0,'alpha',0,'m',m6,'r',[xc6 yc6 0]','I',I6);

% Version from Nicolescu, 2015
% R1 = Link('modified','revolute','alpha',0,'a',0,'d',0.680,'m',m1,'r',[xc1 yc1 zc1]','I',I1);
% R2 = Link('modified','revolute','alpha',-pi/2,'a',0.320,'d',0,'m',m2,'r',[xc2 yc2 zc2]','I',I2);

% Version from my own illustration (see OneNote)
R1 = Link('modified','revolute','alpha',0,'a',0,'d',0.680,'m',m1,'r',[xc1 yc1 zc1]','I',I1);
R2 = Link('modified','revolute','alpha',pi/2,'a',0.320,'d',0,'m',m2,'r',[xc2 yc2 zc2]','I',I2);
R3 = Link('modified','revolute','alpha',pi,'a',0.975,'d',0,'m',m3,'r',[xc3 yc3 zc3]','I',I3);

ABB6620 = SerialLink([R1 R2 R3]);
T30 = ABB6620.A(1:3,[0 30 0]./57.3);
% ABB6620.fkine([0 50 50]./57.3)
% VEL = ABB6620.jacob0([50 50 50 50 50 50]./57.3)*[200 200 200 200 200 200]'./57.3
ABB6620.plot([0 pi/2 0]);

Torq = ABB6620.rne([0 50 50]./57.3,[10 10 10]./57.3,[10 10 10]./57.3);
% Torq = subs(Torq,[m1 xc1 yc1 zc1],[220.489 0.419 0.033 -0.151]);
% Torq = subs(Torq,[Ixy1 Iyz1 Ixz1 Ixx1 Iyy1 Izz1],[1.826 -0.176 -0.543 14.744 25.283 17.320]);
Torq = subs(Torq,[m2 xc2 yc2 zc2],[110.698 0.418 0 0.218]);
Torq = subs(Torq,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[-0.011 -0.005 11.122 6.208 35.317 30.285]);
Torq = subs(Torq,[m3 xc3 yc3 zc3],[250.440 0.156 0.268 -0.001]);
Torq = subs(Torq,[Ixy3 Iyz3 Ixz3 Ixx3 Iyy3 Izz3],[13.428 -1.510 0.559 57.989 11.972 63.456]);
% Torq = subs(Torq,[a1 d1 a2 g],[0.32 0.68 0.975 9.81]);
Torq = vpa(eval(Torq)',5)

% figure;
% hold on;
% for t2 = 0:150
%     Torq = ABB6620.rne([0 t2]./57.3,[0 t2]./57.3,[0 t2]./57.3);
%     Torq = subs(Torq,[m1 xc1 yc1 zc1],[220.489 0.419 0.033 -0.151]);
%     Torq = subs(Torq,[Ixy1 Iyz1 Ixz1 Ixx1 Iyy1 Izz1],[1.826 -0.176 -0.543 14.744 25.283 17.320]);
%     Torq = subs(Torq,[m2 xc2 yc2 zc2],[361.140 0.912 0.186 -0.068]);
%     Torq = subs(Torq,[Ixy2 Iyz2 Ixz2 Ixx2 Iyy2 Izz2],[78.943 -1.515 -10.820 64.197 361.474 407.926]);
%     Torq = eval(Torq);
%     scatter(t2,Torq(2));
% end
% FunTorq = matlabFunction(Torq);