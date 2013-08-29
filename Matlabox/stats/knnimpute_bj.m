function [imputed] = knnimpute_bj(dbase, queries, k, varargin)
    %impute missing values using knn
    %input:
    %   dbase: data to find the knn from, sample x observation (gene)
    %   queries: data to impute, sample x query
    %   k: for knn
    %   median: {true, false}; if false (default), use mean
    %return:
    %   imputed: same size as queries, filled in with imputed values
    %
    
    if nargin < 3,  k = 5;   end
    median = false;
    euclid = false;
    
    i = 1;
    while i < length(varargin)
        switch lower(varargin{i})
            case 'median'
                median = varargin{i+1};
            case 'mean'
                median = ~varargin{i+1};
            case 'euclid'
                euclid = varargin{i+1};
            case 'corr'
                euclid = ~varargin{i+1};
            otherwise
                error('Unknown option %s\n',varargin{i});
        end
        i = i + 2;
    end
    
    [nsamp nobs] = size(dbase);
    [nsamp1 nquery] = size(queries);
    if nsamp ~= nsamp1
        error('inconst sample size\n');
    end
    
    if median
        fhandle = @nanmedian;
    else
        fhandle = @nanmean;
    end
    
    if euclid
        distfun = @eucliddist;
    else
        distfun = @corrdist;
    end
    
    imputed = queries;
    tic;
    for i = 1:nquery
        dist = distfun(dbase, queries(:,i), nobs);
        [tmp si] = sort(dist);        
        clear dist
        nani = isnan(queries(:,i));
        knnavg = fhandle(dbase(:,si(1:k))');        
        imputed(nani,i) = knnavg(nani);
    end
    toc;
    
end

function dist = eucliddist(dbase, query, n)
    dist = nansum((dbase - repmat(query, 1, n)).^2);
end

function dist = corrdist(dbase, query, n)
    dist = - corr(dbase, query, 'rows', 'pairwise');
end