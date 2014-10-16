function histdata = histWithCount( xdata, countdata, binrange, ignoreOutsideBin)
    % xdata and countdata have the same dimension
    % xdata represents the values that will be used to bin the data
    % countdata represents the values that will be used to calculate the
    % histogram
    % eg. xdata = [1, 1, 2, 4]
    %     countdata = [10, 2, 1, 1]
    % if we set the bin to 1:4, the function calculate
    % bin = 1:4, and histdata = [10+2, 1, 0, 1]
    %
    % similarly, accoding to the bin setup, xdata is used to decide which
    % data points in countdata belong to the bin, and those data points are
    % sum up together
    %
    % if countdata are all ones, then this is just as same as hist()
    %
    if nargin < 4, 
        ignoreOutsideBin = false;
    end
    
    if size(xdata, 1) == 1
        xdata = xdata';
    end
    
    %histc assign xdata in [left, right): xdata >= left & xdata < right;
    %except the last bin, which equals to the last number
    [~, idx] = histc(xdata, binrange);
    if any(idx(:) == 0) && ~ignoreOutsideBin
        warning('some data are not in binrage');
    end
    
    [m, n] = size(xdata);
    binsize = length(binrange);
    histdata = zeros(binsize, n);
    %idx is the same size as xdata and countdata
    %idx is now the index to histdata's first dimention (which bin the data
    %goes to)
    dataidx = repmat(1:n, m, 1);
    idx = idx(:);
    dataidx = dataidx(:);
    countdata = countdata(:);
    valid = idx ~= 0;
    histdata(1:max(idx(:)),:) = accumarray( [idx(valid), dataidx(valid)], countdata(valid)); %give max(idx(:)) x n
    