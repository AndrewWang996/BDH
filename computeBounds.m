% function [lowerBoundOnEachSegment, upperBoundOnEachSegment, lowerBoundGlobal, upperBoundGlobal, minOnSamples, maxOnSamples] = computeBounds(L_segments, samples, f)
% %computeBounds - Compute an upper and lower bounds for a real function f on a polygon.
% %
% %   L_segments - the Lipschitz constant of the function f on each segment.
% %   samples - position of the vertices of the polygon. Each two consecutive vertices are considered as one segment (edge) of the polygon.
% %             the first segment endpoints are samples(1) and samples(2).
% %   f - the known values of the function (must be real) at the samples (vertices of the polygon).
% 
% 
% 
%     h = abs(circshift(samples, -1) - samples)/2; %fill distance for each segment
%     
%     avg_val = 1/2*(f + circshift(f, -1));
%     deviation = L_segments.*h;
%     
%     lowerBoundOnEachSegment = avg_val - deviation;
%     upperBoundOnEachSegment = avg_val + deviation;
% 
%     lowerBoundGlobal = min(lowerBoundOnEachSegment);
%     upperBoundGlobal = max(upperBoundOnEachSegment);
% 
%     minOnSamples = min(f);
%     maxOnSamples = max(f);
%     
% end




