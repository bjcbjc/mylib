function fns = listfilename(pattern, usels)
    if nargin < 2, usels = false; end
    if ~usels
        f = dir(pattern);
        fns = cell(size(f));
        for i = 1:length(fns)
            fns{i} = f(i).name;
        end
    else
        f = ls(pattern);
        fns = textscan(f, '%s', 'delimiter', ' \t\n');
        fns = fns{1};
    end
    fns(cellfun(@isempty, fns)) = [];
    fns = sort(fns);