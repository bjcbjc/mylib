function [dataidx, thresholds] = getOutlier(data, method, varargin)

    para.percentileRange = [2.5, 97.5];
    para.pCut = 0.001;

    para = assignpara(para, varargin{:});
    
    if ~isvector(data)
        error('data has to be a vector');
    end
    
    thresholds = NaN(1,2);
    dataidx = [];
    if strcmpi(method, 'expDiffEstimate')
%         q1 = prctile(data, para.percentileRange(1));
%         q2 = prctile(data, para.percentileRange(2));
%         valididx = data >= q1 & data <= q2;
        
        %estimate parameter for exponential distribution, which is used to
        %estimate the distribution of difference between sorted data values
        sdata = sort(data);
        g = gradient(sort(data));
        q1 = prctile(g, para.percentileRange(1));
        q2 = prctile(g, para.percentileRange(2));
        mu = mean(g( g >= q1 & g <= q2));
        %mu = mean(gradient(sort(data(valididx))));
        pval = 1 - expcdf( gradient( sdata ), mu);
        m = mean(data);
        
        outliers = sdata(pval < para.pCut);
        minv = max(outliers(outliers < m));
        maxv = min(outliers(outliers > m));
        dataidx = false(length(data), 1);
        if ~isempty(minv)
            dataidx = dataidx | data < minv;
            thresholds(1) = minv;
        end
        if ~isempty(maxv)
            dataidx = dataidx | data > maxv;
            thresholds(2) = maxv;
        end
        dataidx = find(dataidx);
    else
        fprintf('unknown method\n');
    end
