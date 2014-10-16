function [varargout] = binCount(data, bin, doplot)
    % left open but right close bin count (different from histc)
    % eg. bin = [0, 5, 10]
    %  count: [  -Inf<x<=0, 0<x<=5, 5<x<=10, 10<x<=Inf ]
    % return length(bin)+1 count
    %  binLabel: {'0', '5', '10', '10+}
    %
    %
    
    if nargin < 3, doplot = false; end
    
    if any(isinf(bin))
        error('bin cannot contain Inf');
    end
    
    bin = sort(bin(:));
    nBin = length(bin) + 1;    
    bin = [-Inf; bin; Inf];
    [m, n] = size(data);
    
    if m > 1 && n > 1
        count = zeros(nBin, n);
    else        
        count = zeros(nBin, 1);               
    end
    
    for i = 1:nBin
        count(i,:) = sum( data > bin(i) & data <= bin(i+1) );
    end
    
    binLabel = numarray2strarray(bin(2:end)); 
    binLabel{end} = [binLabel{end-1}, '+'];
    
    if doplot
        handle = bar(count);
        set(gca, 'xtick', 1:nBin, 'xticklabel', binLabel);
        xlim([0 nBin+1]);
    else
        handle = 0;
    end
    
    output = {count, binLabel, handle};
    for i = 1:nargout
        varargout{i} = output{i};
    end
    