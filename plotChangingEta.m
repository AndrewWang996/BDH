%{


%}
function w = plotChangingEta(eta, nsteps)

n = length(eta);
t = linspace(0,1,n);
tt = linspace(0,1,(n-1)*nsteps + 1);
y = [transpose(real(eta)) ; transpose(imag(eta))];
etaI = spline(t, y);
yy = ppval(etaI,tt);
plot(y(1,:),y(2,:),'o',yy(1,:),yy(2,:),'o'), axis equal

end
