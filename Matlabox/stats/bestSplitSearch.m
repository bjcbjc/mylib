function [bestscore, splitvalue, splittype] = bestSplitSearch(Y, X, sizelimit, testfun)
    % search for best binary split point
    %
    % Y: #samples x #variables
    % X: #samples x #variables
    % sizelimit: minimum number of samples in one side
    % testfun: function handle
    %
    
    if nargin < 4
        testfun = @normalgamma;
    end
    
    if nargin < 3
        sizelimit = 5;
    end
    
    [n1, nx] = size(X);
    [n2, ny] = size(Y);
    
    if n1 ~= n2
        error('X and Y must have the same number of samples');
    end
    
    bestscore = NaN(ny, nx);
    splitvalue = NaN(ny, nx);
    splittype = - ones(ny, nx); %+ for >= and - for <=
    
    sortX = sort(X, 1);
    %cut the range based on sizelimit
    sortX(1:(sizelimit-1), :) = [];
    sortX(end-(sizelimit-2):end, :) = [];
    
    for xi = 1:nx
        uvalue = unique(sortX(:, xi));
        nsplit = length(uvalue);
        testscore = -Inf(ny, nsplit);
        for splitidx = 1:nsplit
            if splitidx > (nsplit / 2)
                sampidx = X(:, xi) >= uvalue(splitidx);
            else
                sampidx = X(:, xi) <= uvalue(splitidx);
            end
            if sum(sampidx) < sizelimit || sum(~sampidx) < sizelimit
                continue;
            end
            minNumValidY = min( sum(~isnan(Y(sampidx, :)), 1), sum(~isnan(Y(~sampidx,:)), 1) );
            validY = minNumValidY >= sizelimit;
            testscore(validY, splitidx) = testfun( Y(sampidx, validY), Y(~sampidx, validY) );
        end
        [bestscore(:, xi), bestidx] = max(testscore, [], 2);
        splitvalue(:, xi) = uvalue(bestidx);
        splittype( bestidx > nsplit/2, xi) = 1;
    end
    
end

function score = normalgamma(data1, data2, alpha, lambda)
    if nargin < 3, alpha = 2; end
    if nargin < 4, lambda = 1; end
    data1 = data1';
    data2 = data2';
    score = NormalGammaScore_ForMatrix(alpha, lambda, data1) + ...
        NormalGammaScore_ForMatrix(alpha, lambda, data2) - ...
        NormalGammaScore_ForMatrix(alpha, lambda, [data1 data2]);
end