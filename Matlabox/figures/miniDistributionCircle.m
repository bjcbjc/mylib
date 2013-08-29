function surfhandle = miniDistributionCircle(centerX, centerY, radius, data)


    data = sort(data(:))';
    n = length(data);
    
    if n < 10
        scale = ceil(10/n) + 1;
        data = sort(repmat(data, 1, scale));
        n = length(data);
    end
    
    [x, y] = pol2cart(linspace(0, 2*pi, n+1), ones(1,n+1) * radius);
    
    x = [zeros(1,n+1); x] + centerX;
    y = [zeros(1,n+1); y] + centerY;
    color = repmat(data, 2, 1); %last row is not used
    color(:, end+1) = color(:, end); %last column is not used
    
    surfhandle = pcolor(x, y, color);
    
