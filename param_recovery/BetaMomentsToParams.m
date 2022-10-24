
function [a, b] = BetaMomentsToParams(m, s)
% Equivalent parameters for a Beta distribution with mean m and sd s
% Apply formula
v = s.^2;
a = ((1-m)./v -1./m).*(m.^2);
b = a.*(1./m - 1);

% round up to a given precision to avoid degeneracy
a = round(a*1e5)/1e5;
b = round(b*1e5)/1e5;
end

