function noise = FunNoise(t,mode)
% return noise needed at each time
% input:
% t - time [s]
% mode - an integer in [1,3] representing noise for position, velocity or
% acceleration
persistent matNoise;
if isempty(matNoise)
    matNoise = randn(1e6,3)/57.3*0.5; % 0.5 deg sigma for position
    matNoise(:,2) = matNoise(:,2)*2; % 1 deg/second sigma for velocity
    matNoise(:,3) = matNoise(:,3)*2; % 1 deg/second^2 sigma for acceleration
end
t_step = 1e-3;
tmp = t/t_step;
k = floor(tmp)+1;
disp(['t = ',num2str(t),' (s)']);
if tmp ~= floor(tmp)
   grad = 1/t_step*(matNoise(k+1,mode)-matNoise(k,mode));
   dt = (tmp-floor(tmp))*t_step;
   noise = matNoise(k,mode) + grad*dt;
elseif tmp == floor(tmp)
    noise = matNoise(k,mode);
end
end