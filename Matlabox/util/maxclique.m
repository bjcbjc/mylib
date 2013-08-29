function bic = maxclique(E, varargin)

    %E: binary, #feature (V) x #cond (V)
    
    para.minv = 2;
    
    para = assignpara(para, varargin{:});
    
    [nu nv] = size(E);
    
    if nu ~= nv
        fprintf('E should be symmetric\n');
        bic = [];
        return
    end
    
    dv = sum(E);
    
    %filter matrix
    E(dv < para.minv, dv < para.minv) = 0;
    
    dv = sum(E, 2);
    
    [tmpc tmpv] = eleCounts(dv);
    [tmpv si] = sort(tmpv, 'descend');
    tmpc = cumsum(tmpc(si));
    maxnv = tmpv(find(tmpc >= para.minv, 1, 'first'));

    %initiate
    bic = false(nv, 1);

    %maxnv:-1:para.minv will be the size of sets to enumerate
    %maxnv guarantees at least there are enough features (>=para.minv) that
    %are potenentially all used in maxnv models
    fprintf('looping from %d to %d v nodes\n', maxnv, para.minv);
    for k = maxnv:-1:para.minv
        %goal: find bicliques that have k from V sets
        if mod(k,5) == 0, fprintf('k = %d ...\n',k); end
        %limit on rows that have at least k edges

        vi = find( sum(E, 2) >= k );
        n = length(vi);
        if n < k, continue; end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %enumerate sets of k size from  vi
        vii = 1:k;
        while 1
            if vii(k) > n
                newvii = vii;
                last = 1;
                while 1
                    j = k - last;
                    if j < 1, break; end
                    if newvii(j) + 1 <= n-last
                        newvii(j:end) = newvii(j)+1 : newvii(j)+(k-j+1);
                        break;
                    else
                        last = last + 1;
                    end
                end
                if all(newvii == vii)
                    break
                else
                    vii = newvii;
                end
            end
            %test
            if ~any( E(vi(vii), vi(vii)) == 0)
                bic(vi(vii)) = true;
                return
            end
            vii(k) = vii(k) + 1;
        end

    end
end

