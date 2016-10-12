function w = splineWeight(eta, t)

len = size(eta, 1);
n = size(eta, 2);
times = linspace(0, 1, n);
w = zeros(len, 1);

for i=1:len
   spl = spline(times, eta(i,:));
   w(i) = ppval(spl, t);
end


end