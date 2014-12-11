function [mtx, concord, total] = pairwiseConcord(data, varargin)
    
    para.excludeMissing = 'all'; %{'all', 'any', 'none'};
    para.upper = 'bin'; %{'bin', 'noref', 'symmetric'};
    para = assignpara(para, varargin{:});
          
%     valid = true(size(data));
%     if para.excludeMissing && nnz(data < 0) > 0
%         fprintf('pairwiseConcord: exclude missing genotype\n');
%         valid = data >= 0;
%         idx = any(data<0, 2);
%         fprintf('pairwiseConcord: %d -> %d\n', length(idx), length(idx)-sum(idx));
%         data(idx, :) = [];
%     end
    
    n = size(data, 2); 
    concord = NaN(n, n);
    total = NaN(n, n);

    upperMode = 2;
    if strcmpi(para.upper, 'bin')
        bindata = data ~= 0;
        upperMode = 0;
    else
        bindata = [];
        if strcmpi(para.upper, 'noref')
            upperMode = 1;            
        end
    end
    
    
    for i = 1:n-1
        for j = i+1:n
            if strcmpi(para.excludeMissing, 'all')
                vi = ~all(data(:, [i, j]) == -1, 2);
            elseif strcmpi(para.excludeMissing, 'any')
                vi = ~any(data(:, [i, j]) == -1, 2);
            else
                vi = true(size(data, 1),1);
            end
            
            %lower: genotype
            concord(j, i) = sum(data(vi,i) == data(vi,j));
            total(j, i) = sum(vi);
            
            if upperMode == 0                
                %upper: variant or not
                concord(i, j) = sum(bindata(vi,i) & bindata(vi,j));
                total(i, j) = sum(bindata(vi,i) | bindata(vi,j));
            elseif upperMode == 1
                %upper: genotype, exclude 0/0 == 0/0 and 0/0 ~= ./.
                vi = vi & any(data(:, [i,j]) > 0, 2);
                concord(i, j) = sum(data(vi,i) == data(vi,j));
                total(i, j) = sum(vi);
            else
                concord(i, j) = concord(j, i);
                total(i, j) = total(j, i);
            end
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