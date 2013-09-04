function n = str2double_fast(c,formatstr)
    if nargin < 2
        formatstr = '%g#';
    end
    
    if ~strcmp(formatstr(end), '#')
        formatstr = strcat(formatstr, '#');
    end
    %for non-complex numbers
    n = sscanf(sprintf('%s#', c{:}), formatstr);
    if numel(n) < numel(c)
        n = str2double(c);
    else
        nnum = length(strfind(formatstr, '%'));
        if nnum > 1 && any(size(c)==1)
            ndim = max(size(c));
            n = reshape(n, nnum, ndim)';
        else
            n = reshape(n, size(c));
        end
    end
    
    
    