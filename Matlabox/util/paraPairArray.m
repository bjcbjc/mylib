function d = paraPairArray(para)

    fds = fieldnames(para);
    nfds = length(fds);
    
    d = cell( nfds*2, 1);
    for i = 1:nfds
        d{ (i-1)*2 + 1 } = fds{i};
        d{ i*2 } = para.(fds{i});
    end