function d = distance2polygon(x, y, X, Y)
fprintf('\n Running distance2polygon');
% x = gx(:,1);
% y = gx(:,2);
% X = px(:,1);
% Y = px(:,2);

m = numel(X); 
n = numel(x);

X = repmat(X(:), 1, n); 
Y = repmat(Y(:), 1, n); 

x = repmat(x(:)', m, 1); 
y = repmat(y(:)', m, 1); 

% projection 
uX = X([2:end 1], :) - X; 
uY = Y([2:end 1], :) - Y; 

LAMBDA = ((x - X) .* uX + (y - Y) .* uY) ./ (uX.^2 + uY.^2);

PX = X + LAMBDA .* uX; 
PY = Y + LAMBDA .* uY;

% distance to projections points which is inside a segment 
PD = NaN(size(x)); 
b = ((LAMBDA >= 0) & (LAMBDA <= 1)); 
PD(b) = sqrt((PX(b) - x(b)).^2 + (PY(b) - y(b)).^2);

% distance to polygon point 
D = sqrt((X - x).^2 + (Y - y).^2); 

% distance 
d = min([D; PD]); 

d = d';

% change sign if the point is inside the polygon 
b = inpolygon(x(1, :), y(1, :), X(:, 1), Y(:, 1)); 
d(b) = -d(b);

% figuredocked;
% fDrawPoly = @(x) drawmesh([1:size(x,1); 2:size(x,1) 1; 1:size(x,1)]', x);
% h = fDrawPoly(px);
% hold on;
% bbox = minmax( [X Y]' );
% i = d < -min(range(bbox,2))/100;
% % i = b;
% hp = drawPointSet(gx(i,:), 2, 'b');
% hp2 = drawPointSet(gx(~i,:), 1, 'r');
