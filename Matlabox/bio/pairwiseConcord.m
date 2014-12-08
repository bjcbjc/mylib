function [mtx, concord, total] = pairwiseConcord(data, varargin)
    
    para.excludeMissing = true;
    para = assignpara(para, varargin{:});
        
    valid = true(size(data));
    if para.excludeMissing && nnz(data < 0) > 0
        fprintf('pairwiseConcord: exclude missing genotype\n');
        valid = data >= 0;
%         idx = any(data<0, 2);
%         fprintf('pairwiseConcord: %d -> %d\n', length(idx), length(idx)-sum(idx));
%         data(idx, :) = [];
    end
    
    n = size(data, 2); 
    concord = NaN(n, n);
    total = NaN(n, n);

    bindata = data ~= 0;
    
    for i = 1:n-1
        for j = i+1:n
            vi = all(valid(:, [i, j]), 2);
            %upper: variant or not
            concord(i, j) = sum(bindata(vi,i) & bindata(vi,j));
            total(i, j) = sum(bindata(vi,i) | bindata(vi,j));
            
            %lower: genotype
            concord(j, i) = sum(data(vi,i) == data(vi,j));
            total(j, i) = sum(vi);
        end
    end
    mtx = concord ./ total;
    
%     for i = 1:n-1        
%         j = i+1:n;
%         
%         %upper: variant or not        
%         concord(i, j) = sum(bsxfun(@and, bindata(:,i), bindata(:,j)), 1);
%         total(i, j) = sum(bsxfun(@or, bindata(:,i), bindata(:,j)), 1);           
%         
%         %lower: genotype
%         concord(j, i) = sum(bsxfun(@eq, data(:,i), data(:,j)), 1);
%         total(j, i) = nvar;
%     end
%     mtx = concord ./ total;