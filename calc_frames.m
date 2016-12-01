function [frames] = calc_frames(DerivativeOfCauchyCoordinates, Phi)
    fprintf('\n Running calc_frames');
    fz = DerivativeOfCauchyCoordinates*Phi;
    assert(all(fz ~= 0));
    frames = abs(fz)./fz;

end
