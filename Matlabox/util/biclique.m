function bic = biclique(E, varargin)

    %E: binary, #feature (U) x #cond (V)
    
    para.minv = 5;
    para.minu = 2;
    
    para = assignpara(para, varargin{:});
    
    [nu nv] = size(E);
    
    if nu < nv
        fprintf('transponse E to loop over smaller dimension\n');
        bic = [];
        return
    end
    
    du = sum(E, 2);
    dv = sum(E);
    
    %filter matrix
    E(:, dv < para.minu) = 0;
    E(du < para.minv, :) = 0;

    
    du = sum(E, 2);
    dv = sum(E);
    
    [tmpc tmpu] = eleCounts(du);
    [tmpu si] = sort(tmpu, 'descend');
    tmpc = cumsum(tmpc(si));
    maxnv = tmpu(find(tmpc >= para.minu, 1, 'first'));

    bic = cell(2, 1);
    %initiate
    bic{1} = false(nu, 1);
    bic{2} = false(nv, 1);
    %maxnv:-1:para.minv will be the size of sets to enumerate
    %maxnv guarantees at least there are enough features (>=para.minu) that
    %are potenentially all used in maxnv models
    fprintf('looping from %d to %d v nodes\n', full(maxnv), para.minv);
    for k = maxnv:-1:para.minv
        %goal: find bicliques that have k from V sets
        if mod(k,5) == 0, fprintf('k = %d ...\n',k); end
        %limit on rows that have at least k edges
        %%%%% should loop until the matrix has met the criteria !!!!!
        vi = 1:nv;
        ui = 1:nu;
        while 1
            ui = ui( sum(E(ui, vi), 2) >= k );
            vi = vi( sum(E(ui, vi)) >= para.minu );
            if all(sum(E(ui, vi), 2) >= k) && all(sum(E(ui, vi)) >= para.minu)
                break
            elseif isempty(ui) || isempty(vi)
                break
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        m = length(ui);
        n = length(vi);        
        if n < k, continue; end
        
        %determine if enumerate sets from U, what is the max size
        dv = sum( E(ui, vi) );
        [tmpc tmpu] = eleCounts(dv);
        [tmpu si] = sort(tmpu, 'descend');
        tmpc = cumsum(tmpc(si));
        maxnu = tmpu(find(tmpc >= k, 1, 'first'));
        
        if nsetenumerate(m, maxnu, para.minu) < 2^log2nChoosek_vec(n,k) %nchoosek(n, k)
            %enumerate sets from ui
            for alpha = maxnu:-1:para.minu
                uii = 1:alpha;
                while 1
                    if uii(alpha) > m
                        newuii = uii;
                        last = 1;
                        while 1
                            j = alpha - last;
                            if j < 1, break; end
                            if newuii(j) + 1 <= m-last
                                newuii(j:end) = newuii(j)+1 : newuii(j)+(alpha-j+1);
                                break;
                            else
                                last = last + 1;
                            end
                        end
                        if all(newuii == uii)
                            break
                        else
                            uii = newuii;                            
                        end                    
                    end
                    %test
                    vii = sum( E(ui(uii), vi) ) == alpha;
                    %uii and vii give complete graph, any subgraph
                    %is a clique; but we only want the max clique;
                    %if the max clique has > k nodes, it is already
                    %included; if it is smaller, skip it
                    nsubv = sum(vii);
                    if nsubv == k % we want k v nodes
                        uivec = false(nu, 1); uivec(ui(uii)) = true;
                        vivec = false(nv, 1); vivec(vi(vii)) = true;
                        usubseti = all(bsxfun(@minus, bic{1}, uivec) > -1);
                        vsubseti = all(bsxfun(@minus, bic{2}, vivec) > -1);
                        if ~any(vsubseti & usubseti) %check if it's a subset
                            bic{1} = [bic{1} uivec];
                            bic{2} = [bic{2} vivec];
                        end
                    end
                    uii(alpha) = uii(alpha) + 1;
                end
            end
        else
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
                uii = sum( E(ui, vi(vii)), 2 ) == k;
                if sum(uii) >= para.minu
                    vivec = false(nv, 1); vivec(vi(vii)) = true;
                    vsubseti = all(bsxfun(@minus, bic{2}, vivec) > -1);
                    uivec = false(nu, 1); uivec(ui(uii)) = true;
                    usubseti = all(bsxfun(@minus, bic{1}, uivec) > -1);
                    if ~any(vsubseti & usubseti) %check if it's a subset
                        bic{1} = [bic{1} uivec];
                        bic{2} = [bic{2} vivec];
                    end
                end                
                vii(k) = vii(k) + 1;
            end
        end
    end
    bic{1}(:,1) = [];
    bic{2}(:,1) = [];
end

function nset = nsetenumerate(n, maxk, mink)
    nset = sum(2.^ log2nChoosek_vec(n, maxk:-1:mink));
%     nset = nchoosek(n, maxk);
%     
%     for i = maxk-1:-1:mink
%         nset = nset + nchoosek(n, i);
%     end
end

function val = log2nChoosek_vec(n,k)
% log2nChosek
%  Input:
%   n: A nonnegative integer.
%   k: A nonnegative integer <= n.
%  Output:
%   log base 2 of n choose k, computed to take advantage of the fact that
%   log grows slowly, so that values can be found even for very large
%   n and k.
%
%   vectorize version; if k is a vector, this function is much faster than
%   wrapping it with a for-loop; when k is scalar, it's almost the same
%   speed (can bit tiny slower than) the scalar version
%

    % Check inputs    
    assert(n >= 0 && all(k >= 0), 'n and k must be positive.')
    assert(all(n >= k), 'n must be >= k.')

    % Compute, using a formula that works even for very large n and k.
    alllog = (log2(n:-1:1))';
    maxk = max(k);
    dist = cumsum(alllog(1:maxk) - alllog(n:-1:n-maxk+1));
    val = dist(k);
end
