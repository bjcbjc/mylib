function res = loadStructData(fn)
    res = load(fn);
    fds = fieldnames(res);
    if length(fds) > 1
        error('This is not a structure: %s', fn);
    end    
    res = res.(fds{1});

    
