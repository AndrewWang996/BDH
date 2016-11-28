
%{
This is assuming that the number of keyframes (n) is 3.

%}
function w = quadraticSplineWeight(t)

fprintf('using quadratic spline weight\n');

w = zeros(3, 1);
w(1) = (1-t)*(1-t);
w(2) = 2*(t*(1-t));
w(3) = t*t;

end