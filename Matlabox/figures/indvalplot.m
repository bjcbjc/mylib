function indvalplot(x, data, varargin)
    %treat each column of data as a group, and plot each data plots along
    %Y-axis; each group is distributed on x-axis between [x-0.4 x+0.4]
    %
    
    pointdist = 0.05;
    maxresolution = 100;
    [n ngroup] = size(data);
    
    if length(x) ~= ngroup
        error('length(x) is not the same as #columns in data');
    end
    
    if size(x,1) > size(x,2)
        x = x'; %row vec
    end
    if length(x) == 1
        range = 0.8;
    else
        range = min(diff(x))*0.8;
    end
    
    scatterX = repmat(x, n, 1);
    
    rounddata = data;
    %round the data on y-axis for speed
    if length(unique(data(:))) > maxresolution
        M = max(data(:));
        m = min(data(:));
        bins = m:(M-m)/maxresolution:M;        
        for bi = 1:length(bins)-1
            rounddata(data>=bins(bi) & data<bins(bi+1)) = bins(bi);
        end
    end
    
    for i = 1:ngroup
        [c u] = eleCounts(rounddata(:,i));
        if length(u) == 1, continue; end
        for j = 1:length(u)
            if ( (c(j)-1) * pointdist ) > range 
                scatterX(rounddata(:,i)==u(j), i) = (x(i)-0.5*range):range/(c(j)-1):(x(i)+0.5*range);
            else                
                scatterX(rounddata(:,i)==u(j), i) = (x(i) - 0.5*(c(j)-1)*pointdist): pointdist: (x(i) +0.5*(c(j)-1)*pointdist);
            end
        end
    end
    plot(scatterX, data, varargin{:});
end