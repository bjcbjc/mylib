function [pvalue, nOverlap] = hygetest(data1, data2)
    % data1 and data2 are binary matrix: observation x variable
    % test P( #overlap >= overlap btw variable A in data1 and variable B in
    % data2); if data2 is empty, do pair-wise test among all variables in
    % data1
    % 
    
    if nargin < 2, data2 = data1; end
    if isempty(data2), data2 = data1; end
    
    [m, n1] = size(data1);
    [m2, n2] = size(data2);
    if m ~= m2
        error('size inconst');
    end
            
    nOverlap = NaN(n1, n2);    
    if ~islogical(data1)
        data1 = data1~=0;
    end
    if ~islogical(data2)
        data2 = data2~=0;
    end
    if n1 < n2        
        for i = 1:n1
            nOverlap(i, :) = sum(bsxfun(@and, data1(:,i), data2), 1);
        end
    else    
        for i = 1:n2
            nOverlap(:, i) = sum(bsxfun(@and, data1, data2(i,:)), 1);
        end
    end
    
    nCase1 = sum(data1, 1);
    nCase2 = sum(data2, 1);
    pvalue = hygecdf(nOverlap-1, m, nCase1, nCase2, 'upper');
