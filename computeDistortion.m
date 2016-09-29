function [k, sigma1, sigma2, tau, arap] = computeDistortion(DerivativeOfCauchyCoordinates, Phi, Psy)
%k - conformal distortion
%sigma1 - largest singular value of the Jacobian
%sigma2 - smallest (signed) singular value of the Jacobian - negative value indicates noninjectiveness
%tau - ismetric distortion max(sigma1, 1/sigma2)
%arap - ARAP distortion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_fz = abs(DerivativeOfCauchyCoordinates*Phi);
abs_fzbar = abs(DerivativeOfCauchyCoordinates*Psy);

k = gather(abs_fzbar ./ abs_fz);

if(nargout > 1)
    sigma1 = gather(abs_fz + abs_fzbar);
end

if(nargout > 2)
    sigma2 = gather(abs_fz - abs_fzbar); %note that sigma2 is actually abs(abs_fz - abs_fzbar) but we remove the abs to identify non injectiveness (sigma2 will be negative)
end

if(nargout > 3)
    tau = gather(max(sigma1, 1./sigma2));
end

if(nargout > 4)
    arap = gather((abs_fz - 1).^2 + abs_fzbar.^2); %this expression is only true when the mapping is locally injective - i.e. when (signed) sigma2 is positive
end
end

