function output = FunSat(x)
% Saturation function for normalized variable x
if x < -1
    output = -1;
elseif x < 1
    output = x;
else
    output = 1;
end
end