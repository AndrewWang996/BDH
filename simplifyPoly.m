function [i, msub] = simplifyPoly(x, n, angleThresh)
fprintf('\n Running simplifyPoly');
assert(~isreal(x));

% angleThresh = 0.3;
% n = 20;

i = 1:numel(x);

while numel(i)>n
    y = x(i);
    a = abs( angle( (y([2:end 1]) - y)./(y - y([end 1:end-1])) ) );
    a = min(a, pi-a);

    [~, j] = sort(a);
    if a(j(1)) >= angleThresh
        break;
    end

    j = j(1:max(1, floor((numel(i)-n)/2)));
    k = 2;
    while k<=numel(j)
        if min( abs(j(k) - j(1:k-1)) ) == 1  % do not remove neighboring points simutaneously
            j = j([1:k-1 k+1:end]);
        elseif a(j(k)) > angleThresh
            j = j(1:k-1);
        else
            k = k+1;
        end
    end

%     i = sort(i(setdiff(1:numel(i), j)));
    i = setdiff(i, i(j));

%     delete(h);    h = fDrawPoly( x(i) ); set(h, 'marker', 'x');
end


if nargout>1
    nx = numel(x);
    nnx = numel(i);
    msub = sparse(i, 1:nnx, 1, nx, nnx);

    for j = setdiff(1:nx, i)
        k = find(i<j, 1, 'last');
        if isempty(k); k = nnx; end;

        k_1 = mod(k,nnx)+1;
        w = abs( (x( i(k_1) ) - x(j))/(x( i(k_1) ) - x( i(k) )) );
        msub(j, [k k_1]) = [w 1-w];
    end
end


% fDrawPoly = @(x) plot( real(x([1:end 1])), imag(x([1:end 1])), '-x');
% % fDrawPoly = @(x)drawmesh([1:size(x,1);2:size(x,1),1;1:size(x,1)]',x);
% figuredocked; h = fDrawPoly(x(i)); %h.EdgeColor='r';
% % hold on; h = fDrawPoly(x(i)); set(h, 'marker', 'x');
% % hold on; h = fDrawPoly( msub*x(i) ); set(h, 'marker', 'x');