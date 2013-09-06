function fns = listfilename(pattern)
    f = dir(pattern);
    fns = cell(size(f));
    for i = 1:length(fns)
        fns{i} = f(i).name;
    end