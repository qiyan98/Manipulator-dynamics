clear;
% T_l = [-1e3:50:1e3];
T_l = 800;
% Omega = [-90:1:-5,-5:1e-1:5,5:1:90]/57.3;
Omega = [-90:0.01:-1,1:0.01:90]/57.3;
% Tem = [30:1:50];
Tem = 40;
data = zeros(10,4);
k = 1;
for j=1:length(T_l)
    for jj=1:length(Omega)
        for jjj=1:length(Tem)
            data(k,1) = FunFriction(T_l(j),Omega(jj),Tem(jjj));
            data(k,2:4) = [T_l(j),Omega(jj),Tem(jjj)];
            k = k + 1;
        end
    end
end
save FricData_omega_only.mat data