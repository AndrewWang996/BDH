function w = splineWeight(eta, t)

numVertices = size(eta, 1);
numKeyframes = size(eta, 2);
times = linspace(0, 1, numKeyframes);

w = zeros(numVertices, 1);

persistent splines;
persistent interpolated;
if isempty( splines )
    splines = struct([]);
end
if isempty(interpolated)
    interpolated = zeros(numVertices, 1);
end

for i=1:numVertices
    if i > size(splines,2)
        splines(i).spl = spline(times, eta(i,:));
    end
    w(i) = ppval(splines(i).spl, t);
end


end