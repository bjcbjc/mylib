function res = movingAvg(vector, window)

    dim = size(vector);
    if dim(1) > 1 && dim(2) > 1
        error('input should be a vector');
    end
    vector = vector(:)'; % row
    
    n = max(dim);
    
    if mod(window, 2) == 0
        window = window + 1;
    end
    
    if n > window
        data = NaN(window, n);
        half = (window - 1) / 2;
        for i = 1:half
            data(i, 1:n-i) = vector(i+1:end);            
        end
        data(half+1, :) = vector;
        for i = 1:half
            data(half+i+1, i+1:n) = vector(1:n-i);
        end
    end    
    res = nanmean(data, 1);
end
