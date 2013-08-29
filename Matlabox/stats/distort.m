function dist = distort(data, kidx, centroid, tau_inv)
    %calculate minimum distortion from kmean clustering's result
    %
    %data: #observation x #feature
    %centroid: #k x #feature
    %kidx: #observation x 1, index of cluster for each observation
    %tau_inv: inverse of covariance of data; #feature x #feature    
    %
    %
    
    K = max(kidx);
    
    [n, p] = size(data);
    
    if nargin < 3 || isempty(centroid)
        centroid = NaN(K, p);
        for i = 1:K
            centroid(i, :) = nanmean(data(kidx == i, :), 1);
        end
    end
    
    if nargin < 4 || isempty(tau_inv)
        tau_inv = pinv( cov(data) );
    end        
    
    D = data - centroid( kidx, :);
    
    d = 0;
    for i = 1:n
        d = d + D(i,:) * tau_inv * D(i,:)';
    end
    
    dist = d ./ n ./ p;
    