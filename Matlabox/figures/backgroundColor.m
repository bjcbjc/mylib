function h = backgroundColor(xrange, yrange, nx, ny, color)
    
    x0 = xrange(1);
    y0 = yrange(1);
    
    w = diff(xrange);
    h = diff(yrange);
    
    segXSize = w / nx;
    segYSize = h / ny;
    
    if size(color, 1) < nx * ny
        error('%d * %d segments, but only %d colors\n', ...
            nx, ny, size(color, 1));
    end
    
    h = NaN(nx, ny);
    for yi = 1:ny
        for xi = 1:nx
            x = x0 + segXSize * (xi - 1);
            y = y0 + segYSize * (yi - 1);
            h(xi, yi) = rectangle('position', [x, y, segXSize, segYSize], ...
                'edgecolor', 'none', ...
                'facecolor', color((yi-1)*nx + xi, :));
        end
    end
    uistack(h(:), 'bottom');
end