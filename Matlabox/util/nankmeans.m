function [indexes, score] = nankmeans(x, k, dis, rep, inic)

    if nargin < 5, inic = NaN; end
    if nargin < 4, rep = 50; end
    if nargin < 3, dis = 'euclid'; end
    
    t = cputime();
    
    [n p] = size(x);
    
    maxiter = 400;
    scx = zeromean_univar_normalization(x);

    indexes = NaN(n,1);
    score = Inf;
    for r = 1:rep
        i = 1;
        if isscalar(inic)
            centroids = scx(randsample(n,k), :);
        else
            centroids = inic;
        end
        rindexes = zeros(n,1);
        while i <= maxiter
            i = i + 1;
            d = distance(scx, centroids, k, n, dis);
            [m rindexes] = nanmin(d,[],2);            
            centroids = recenter(rindexes, scx, k, p);
            repscore = sum(m);            
            if abs(repscore-score) <= eps
                break
            end
            if repscore < score
                score = repscore;
                indexes = rindexes;
            end
        end
    end
    fprintf('time: %f\n',cputime()-t);
end

function c = recenter(indexes, x, k, p)
    c = zeros(k, p);
    for i = 1:k
        c(i,:) = nanmean(x(indexes==i,:));
    end
end

function dist = distance(x, c, k, n, method)
    dist = NaN(n, k);
    switch method
        case 'pearson'
            dist = 1-corr(x',c','rows','pairwise');
        case 'spearman'
            dist = 1-corr(x',c','rows','pairwise','type','Spearman');
        otherwise
            for i = 1:k
                if length(find(isnan(c(i,:))==0)) > 1
                    dist(:,i) = nansum(((x - repmat(c(i,:),n,1)).^2)')';
                end
            end
    end
end


