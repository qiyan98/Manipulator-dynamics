function output = signum(s,g)
% Input: state s and parameter g
% Output: sig(s)^g
output = abs(s)^g*sign(s);
end