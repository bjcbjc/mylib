function pred = splitPredict(Y, x, splittype, splitvalue)
    %
    % Y: #sample x #variable
    % X: #sample x 1
    % splittype: #variable x 1; +1 for >= and -1 for <=
    % splitvalue: #variable x 1, thresholds in X
    %
    
    [n, ny] = size(Y);
        
    group = false(n, ny);
    geSplitY = splittype > 0;
    group(:, geSplitY) = ...
        bsxfun(@ge, x, splitvalue(geSplitY)'); %#sample x #geSplitY
    group(:, ~geSplitY) = ...
        bsxfun(@le, x, splitvalue(~geSplitY)'); %#sample x #leSplitY
    
    tmpY = Y;
    tmpY(~group) = NaN;
    pred = repmat(nanmean(tmpY,1), n, 1);
    
    tmpY = Y;
    tmpY(group) = NaN;
    tmpY = repmat(nanmean(tmpY,1), n ,1);
    pred(~group) = tmpY(~group);
    
    
    
    