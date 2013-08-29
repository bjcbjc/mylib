function u = mynanunique(data)

    if ~any(isnan(data(:)))
        u = unique(data);
    else
        u = unique(data(~isnan(data)));
    end