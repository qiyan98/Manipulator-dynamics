function  [t_ctrl,X_ctrl] = discrete_ctrl(EOM,MaxTime,X_0)
% customized ode solver for discrete-time controller implementations

%% controller parameters
dt = 0.005;
size_X = length(X_0); % in case of 3-state PID
t_ctrl = zeros(MaxTime/dt,1);
X_ctrl = zeros(MaxTime/dt,size_X);

%% discrete-time solver
for j = 1:length(t_ctrl)
    cur_t = dt*j;
    if j > 1
        cur_X = X_ctrl(j-1,:);
    else
        cur_X = reshape(X_0,[1 size_X]);
    end
    cur_X_dot = EOM(cur_t,cur_X); 
    cur_X_dot = cur_X_dot'; % to make it row vector
    if j > 1
        X_ctrl(j,:) = cur_X_dot*dt + X_ctrl(j-1,:);
    else
        X_ctrl(j,:) = cur_X_dot*dt + reshape(X_0,[1 size_X]);
    end
    t_ctrl(j) = cur_t;
end


end