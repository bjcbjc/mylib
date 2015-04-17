function [mtx, concord, total] = pairwiseConcord2(data1, data2, varargin)
    %calculate concordance between columns of data1 and columns of data2
    para.excludeMissing = 'all'; %{'all', 'any', 'none'};
    para.excludeRef = false; 
    para = assignpara(para, varargin{:});
          
    n1 = size(data1, 2); 
    n2 = size(data2, 2);
    concord = NaN(n1, n2);
    total = NaN(n1, n2);

    for i = 1:n1
        for j = 1:n2
            if strcmpi(para.excludeMissing, 'all')
                vi = data1(:, i) ~= -1 | data2(:, j) ~= -1;
            elseif strcmpi(para.excludeMissing, 'any')
                vi = data1(:, i) ~= -1 & data2(:, j) ~= -1;
            else
                vi = true(size(data1, 1),1);
            end
            
            if para.excludeRef
                %genotype, exclude 0/0 == 0/0 and 0/0 ~= ./.
                vi = vi & ( data1(:, i) > 0 | data2(:, j) > 0);
            end
            
            %genotype
            concord(i, j) = sum(data1(vi,i) == data2(vi,j));
            total(i, j) = sum(vi);
        end
    end
    mtx = concord ./ total;
    
