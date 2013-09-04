function TF = isnumericstring(cellarrays)

    if ischar(cellarrays)
        cellarrays = {cellarrays};
    end

    [m n] = size(cellarrays);
    TF = false(m, n);
    
    cellarrays = strrep(cellarrays, '.', '');
    cellarrays = strrep(cellarrays, '+', '');
    cellarrays = strrep(cellarrays, '-', '');
    for i = 1:m
        for j = 1:n
            if isempty(cellarrays{i,j})
                TF(i,j) = false;
            elseif ~strcmpi(cellarrays{i,j}, 'nan')
                TF(i,j) = all(cellarrays{i,j} >= '0' & cellarrays{i,j} <= '9');
            else
                TF(i,j) = true;
            end
        end        
    end