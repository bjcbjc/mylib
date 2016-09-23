function fig = offFigure(width, height)
if nargin < 2
    height = width(2);
    width = width(1);
end
fig = figure('position', [0, 0, width, height], 'visible', 'off','paperpositionmode','auto');