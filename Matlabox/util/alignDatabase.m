function newdb = alignDatabase(dbstruct, fd2align, alignlist)
    %align/order data in dbstruct
    %
    % input: 
    %   dbstruct, structure containing the database
    %   fd2align: string, the fieldname to align
    %   alignlist: cell array, list to align to
    %
    % Note: for those that are not in alignlist, create NaN in the matrix
    %
    %
    
    

    [~, i] = ismember(alignlist, dbstruct.(fd2align));
    
    if all(i==0)
        error('nothing matched! Use your brain!');
    end
    
    
    fds = fieldnames(dbstruct);
    nfd = length(fds);
    dim2change = zeros(nfd,1);
    ndim = zeros(nfd,1);
    
    matchdim = length(dbstruct.(fd2align));
    for j = 1:nfd
        dim = size(dbstruct.(fds{j}));        
        ndim(j) = length(dim);
        for di = 1:ndim(j)
            if dim(di) == matchdim                
                dim2change(j) = di;
                break
            end
        end
    end
    
    
    newlen = length(alignlist);
    for j = 1:nfd
        if dim2change(j) == 0
            newdb.(fds{j}) = dbstruct.(fds{j});
        else
            if strcmp(fds{j}, fd2align)
                newdb.(fds{j}) = alignlist;
            else
                ifsparse = issparse(dbstruct.(fds{j}));
                permorder = [dim2change(j) setdiff(1:ndim(j), dim2change(j))];
                backorder = NaN(1, ndim(j));
                backorder(dim2change(j)) = 1;
                backorder(setdiff(1:ndim(j), dim2change(j))) = 2:ndim(j);
                d = permute(dbstruct.(fds{j}), permorder);
                dim = size(dbstruct.(fds{j}));
                newdim = dim(permorder);
                newdim(1) = newlen;
                if iscell(d)
                    newdb.(fds{j}) = cell(newdim);
                else
                    newdb.(fds{j}) = NaN(newdim);
                end
                s = '';
                for k = 1:ndim-1
                    s = [s ',:'];
                end
                s = [s ')'];
                eval(sprintf('newdb.(fds{j})(i~=0%s = d(i(i~=0)%s;',s, s));
                newdb.(fds{j}) = permute(newdb.(fds{j}), backorder);
                if ifsparse
                    newdb.(fds{j}) = sparse(newdb.(fds{j}));
                end
            end
        end
    end
    