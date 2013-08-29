function [res bound] = discdata(data, bin, varargin)
    %discretize data
    %
    %data: #sample x #variable (eg, gene)
    %bin: # of bins
    %varargin:
    %   'method': {'bucket','std','soft'}
    %   'soft': factor for +-std on the boundary; def = 0.2;
    %   'bound': use it to discretize the data
    %
    
    if nargin < 2, bin =  5; end
    
    para.method = 'std';    
    para.soft = 0.2;
    para.bound = [];
    para.stdfactor = [];
    para = assignpara(para, varargin{:});
    
    [para.nsamp para.nvar] = size(data);
    
    if isempty(para.stdfactor)
        para.stdfactor = setdiff(-(bin-1)/2:(bin-1)/2,0)';
    elseif size(para.stdfactor,1) ~= bin-1 && ...
            (strcmpi(para.method, 'std') || strcmpi(para.method, 'soft'))
        error('length(stdfactor) ~= #bin-1');
    end
    
    %stdfactor should be either a column (#bin-1 x 1) or a matrix (#bin-1 x #gene)
    if size(para.stdfactor,1) == 1, para.stdfactor = para.stdfactor'; end 
    
    if ~isempty(para.bound)
        res = disc(data, bin, para);
        bound = para.bound;
    elseif strcmpi(para.method, 'bucket')
        [res bound] = bucket(data, bin, para);
    elseif strcmpi(para.method, 'softbucket')
        [res lbound ubound] = softbucket(data, bin, para);
        bound{1} = lbound;
        bound{2} = ubound;
    elseif strcmpi(para.method, 'std')
        [res bound] = stddisc(data, bin, para);
    elseif strcmpi(para.method, 'soft')
        [res lbound ubound] = softdisc(data, bin, para);
        bound{1} = lbound;
        bound{2} = ubound;
    else
        error('Unknown method %s\n',para.method);
    end    
    %double check
    tmp = unique(res(:));
    tmp = tmp(~isnan(tmp));
    if nnz(tmp<1) || nnz(tmp>bin)
        error('something is wrong, invalid disc values');
    end
end

function res = disc(data, bin, para)
    res = NaN(para.nsamp, para.nvar);
    if strcmpi(para.method, 'std') || strcmpi(para.method, 'bucket')
        for i = 1:para.nvar
            for bi = 1:bin
                index = data(:,i) > para.bound(bi,i) & data(:,i) <= para.bound(bi+1,i);
                res(index,i) = bi;
            end
        end
    else
        for i = 1:para.nvar
            for bi = 1:bin
                index = data(:,i) > para.bound{1}(bi,i) & data(:,i) < para.bound{2}(bi,i);
                res(index,i) = bi;
            end
            for bi = 1:bin-1
                index = data(:,i) >= para.bound{2}(bi,i) & data(:,i) <= para.bound{1}(bi+1,i);
                res(index,i) = ( data(index,i) - para.bound{2}(bi,i) ) / ...
                    ( para.bound{1}(bi+1,i) - para.bound{2}(bi,i) ) + bi;
            end
        end
    end
end

function [res bound] = bucket(data, bin, para)
    [tmp si] = sort(data); %NaN is sorted at the end
    res = NaN(para.nsamp, para.nvar);
    bound = [-Inf*ones(1, para.nvar); NaN(bin-1, para.nvar); Inf*ones(1, para.nvar)];
    for i = 1:para.nvar
        c = 0;
        bsize = floor(sum(~isnan(data(:,i)))/bin);
        bres = mod(sum(~isnan(data(:,i))), bin);
        for bi = 1:bin
            if bi <= bres
                res(si(c+1:c+bsize+1,i),i) = bi;
                if bi < bin, bound(bi+1,i) = tmp(c+bsize+1,i); end
                c = c + bsize + 1;
            else
                if c+bsize == 0, continue; end
                res(si(c+1:c+bsize,i),i) = bi;
                if bi < bin, bound(bi+1,i) = tmp(c+bsize,i); end
                c = c + bsize;
            end
        end
    end    
end

function [res lbound ubound] = softbucket(data, bin, para)
    [tmp si] = sort(data);    
    
    res = NaN(para.nsamp, para.nvar);
    lbound = [-Inf * ones(1, para.nvar); NaN(bin-1, para.nvar)]; %#bin x nvar
    ubound = [NaN(bin-1, para.nvar); Inf * ones(1, para.nvar)];
    for i = 1:para.nvar
        %calculate bin size and size of soft-region
        gnsamp = sum(~isnan(data(:,i)));
        bsize = [floor(gnsamp*(1-para.soft)/bin), ...
            repmat(floor(gnsamp*(1-2*para.soft)/bin), 1, bin-2), ...
            floor(gnsamp*(1-para.soft)/bin)];
        softtotal = gnsamp-sum(bsize);
        nbiggersoft = mod(softtotal, bin-1);
        softsize = [repmat(ceil(softtotal/(bin-1)), 1, nbiggersoft), ...
            repmat(floor(softtotal/(bin-1)), 1, bin-1-nbiggersoft) ];

        countindex = zeros(1, 2*bin);
        countindex(2:2:2*bin) = bsize;
        countindex(3:2:2*bin-1) = softsize;
        countindex = cumsum(countindex);
        countindex(1:2:2*bin-1) = countindex(1:2:2*bin-1) + 1;
        %now countindex will be [a1, b1, a2, b2...]
        %each pair of (a, b) is the indexes of the sorted element for each bin
        %eg, sorted_data(a1):sorted_data(b1) are the range of bin#1
        for bi = 1:bin
            res( si( countindex(bi*2-1) : countindex(bi*2) ) , i) = bi;
            if bi ~= bin
                ubound(bi, i) = tmp(countindex(bi*2), i);
            end
            if bi ~= 1
                lbound(bi, i) = tmp(countindex(bi*2-1), i);
            end
        end
        for bi = 1:bin-1 %for soft-region
            res( si( countindex(bi*2)+1: countindex(bi*2+1)-1 ), i) = ...
                (tmp( countindex(bi*2)+1: countindex(bi*2+1)-1 ) - ubound(bi) ) / ...
                (lbound(bi+1) - ubound(bi)) + bi;
        end
    end    
end

function [res bound] = stddisc(data, bin, para)
    if mod(bin,2) == 0
        error('#bin should be odd')
    end
    s = nanstd(data);
    m = nanmean(data);
    res = NaN(para.nsamp, para.nvar);
    if size(para.stdfactor, 2) == 1
        bound = [-Inf * ones(1, para.nvar); ...
            repmat( para.stdfactor, 1, para.nvar ); ...
            Inf * ones(1, para.nvar)]; %bin+1 x nvar        
    else
        bound = [-Inf * ones(1, para.nvar); ...
            para.stdfactor; ...
            Inf * ones(1, para.nvar)]; %bin+1 x nvar
    end
    bound = repmat(m, bin+1, 1) + bound .* repmat(s, bin+1, 1);
    for i = 1:para.nvar
        for bi = 1:bin
            index = data(:,i) > bound(bi,i) & data(:,i) <= bound(bi+1,i); %consistent with bucket
            res(index,i) = bi;
        end
    end
end

function [res lbound ubound] = softdisc(data, bin, para)
    if mod(bin,2) == 0
        error('#bin should be odd')
    end
    s = nanstd(data);
    m = nanmean(data);
    res = NaN(para.nsamp, para.nvar);
    if size(para.stdfactor, 2) == 1
        lbound = [-Inf * ones(1, para.nvar); ...
            repmat( para.stdfactor + para.soft, 1, para.nvar )];
        ubound = [repmat( para.stdfactor - para.soft, 1, para.nvar ); ...
            Inf * ones(1, para.nvar)]; %bin x nvar
    else
        lbound = [-Inf * ones(1, para.nvar); ...
            para.stdfactor + para.soft ];
        ubound = [ para.stdfactor - para.soft; ...
            Inf * ones(1, para.nvar)]; %bin x nvar
    end
    lbound = repmat(m, bin, 1) + lbound .* repmat(s, bin, 1);
    ubound = repmat(m, bin, 1) + ubound .* repmat(s, bin, 1);
    for i = 1:para.nvar
        for bi = 1:bin
            index = data(:,i) > lbound(bi,i) & data(:,i) < ubound(bi,i);
            res(index,i) = bi;            
        end
        for bi = 1:bin-1
            index = data(:,i) >= ubound(bi,i) & data(:,i) <= lbound(bi+1,i);
            res(index,i) = (data(index,i)-ubound(bi,i))/(lbound(bi+1,i)-ubound(bi,i)) + bi;
        end
    end
end