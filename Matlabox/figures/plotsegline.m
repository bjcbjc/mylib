function plotsegline(ticks, xyseg, varargin)
    %ticks: position to plot segment lines
    %xyseg: [true/false, true/false], plot x (vertical) and/or y (horizontal) segment lines
    %varargin: pass to plot
    %

    if nargin < 2, xyseg = [true, false]; end
    
    %make ticks to column vec
    ticks = ticks(:);
    
    x = [];
    y = [];
        
    if xyseg(1) %plot vertical segment lines
        x = [x repmat(ticks', 2, 1)];
        y = [y repmat(ylim', 1, length(ticks))];    
    end
    if xyseg(2) %plot horizontal segment lines
        x = [x repmat(xlim', 1, length(ticks))];
        y = [y repmat(ticks', 2, 1)];
    end
    
    hold on;
    line(x, y, varargin{:});
    hold off;