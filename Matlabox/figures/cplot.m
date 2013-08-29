function cplot(x1, x2, varargin)

[~, density, x, y] = kde2d([x1 x2], 64);
%figure
contour(x, y, density, 12, varargin{:});