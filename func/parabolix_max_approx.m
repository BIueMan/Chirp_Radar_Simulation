% find the maximum from a parabola that is build from 3 point
function [x_vertex, y_vertex] = parabolix_max_approx(x, y, warn)
    if nargin < 3
        warn = true;  % Display warning by default
    end
    % corect input
    if size(x, 2) ~= 3 || size(y, 2) ~= 3
        error('Input matrices must be of size m x 3.')
    end
    if warn && (size(x, 1) == 3 || size(y, 1) == 3)
        warning('Input matrices is 3x3, make sure they are not fliped')
    end
    
    function [A, B, C] = parabola_vertex(x, y)
        denom = (x(:,1)-x(:,2)) .* (x(:,1)-x(:,3)) .* (x(:,2)-x(:,3));

        A     = (x(:,3) .* (y(:,2)-y(:,1)) + x(:,2) .* (y(:,1)-y(:,3)) + x(:,1) .* (y(:,3)-y(:,2))) ./ denom;
        B     = (x(:,3).*x(:,3) .* (y(:,1)-y(:,2)) + x(:,2).*x(:,2) .* (y(:,3)-y(:,1)) + x(:,1).*x(:,1) .* (y(:,2)-y(:,3))) ./ denom;
        C     = (x(:,2) .* x(:,3) .* (x(:,2)-x(:,3)) .* y(:,1)+x(:,3) .* x(:,1) .* (x(:,3)-x(:,1)) .* y(:,2)+x(:,1) .* x(:,2) .* (x(:,1)-x(:,2)) .* y(:,3)) ./ denom;
    end

    [a, b, c] = parabola_vertex(x, y);
    
    x_vertex = -b./(2*a);
    y_vertex = a.*x_vertex.^2 + b.*x_vertex + c;
end