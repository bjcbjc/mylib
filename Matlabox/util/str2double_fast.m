function n = str2double_fast(c)
    %for non-complex numbers
    n = sscanf(sprintf('%s#', c{:}), '%g#');
    if numel(n) < numel(c)
        n = str2double(c);
    else
        n = reshape(n, size(c));
    end